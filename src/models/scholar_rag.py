import google.generativeai as genai
from google.cloud import aiplatform
from vertexai.language_models import TextEmbeddingModel
from langchain_community.vectorstores import FAISS
from langchain_core.embeddings import Embeddings
from typing import Dict, Any, List
import logging
import os
import traceback
from dotenv import load_dotenv
from ..services.pubmed_service import PubMedService

# Initialize logging with more details
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Load environment variables
load_dotenv()

# Configure Google credentials
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")
GOOGLE_CLOUD_PROJECT = os.getenv("GOOGLE_CLOUD_PROJECT", "medical-search-2025")

if not GOOGLE_API_KEY:
    raise ValueError("GOOGLE_API_KEY not found in environment variables")

# Configure Gemini
try:
    logger.info(f"Configuring Gemini with project: {GOOGLE_CLOUD_PROJECT}")
    genai.configure(api_key=GOOGLE_API_KEY, project=GOOGLE_CLOUD_PROJECT)
    
    # List available models as a diagnostic step
    models = genai.list_models()
    available_models = [model.name for model in models]
    logger.info(f"Available models: {available_models}")
    
    if "models/gemini-2.0-flash-001" not in available_models:
        logger.warning("models/gemini-2.0-flash-001 not found in available models!")
except Exception as e:
    logger.error(f"Error configuring Gemini: {e}")
    traceback.print_exc()

# Initialize Google Cloud SDK
try:
    aiplatform.init(project=GOOGLE_CLOUD_PROJECT, location="us-central1")
    logger.info(f"Initialized Google Cloud SDK with project: {GOOGLE_CLOUD_PROJECT}")
except Exception as e:
    logger.warning(f"Could not initialize Google Cloud SDK: {e}")
    traceback.print_exc()

class CustomVertexEmbeddings(Embeddings):
    """Custom wrapper for Vertex AI Embeddings to avoid compatibility issues."""
    
    def __init__(self, model_name: str = "text-embedding-005"):
        """Initialize with the model_name."""
        logger.info(f"Initializing CustomVertexEmbeddings with model: {model_name}")
        try:
            self.model = TextEmbeddingModel.from_pretrained(model_name)
            logger.info("TextEmbeddingModel initialization successful")
        except Exception as e:
            logger.error(f"Error initializing TextEmbeddingModel: {e}")
            traceback.print_exc()
            raise
    
    def embed_documents(self, texts: List[str]) -> List[List[float]]:
        """Embed documents using Vertex AI."""
        try:
            embeddings = self.model.get_embeddings(texts)
            logger.debug(f"Generated embeddings for {len(texts)} documents")
            return [embedding.values for embedding in embeddings]
        except Exception as e:
            logger.error(f"Error embedding documents: {e}")
            traceback.print_exc()
            raise
    
    def embed_query(self, text: str) -> List[float]:
        """Embed a query using Vertex AI."""
        try:
            embeddings = self.model.get_embeddings([text])
            logger.debug("Generated embedding for query")
            return embeddings[0].values
        except Exception as e:
            logger.error(f"Error embedding query: {e}")
            traceback.print_exc()
            raise

class ScholarRAG:
    def __init__(self):
        logger.info("Initializing ScholarRAG")
        try:
            self.model = genai.GenerativeModel('gemini-2.0-flash-001')
            logger.info("GenerativeModel initialization successful")
        except Exception as e:
            logger.error(f"Error initializing GenerativeModel: {e}")
            traceback.print_exc()
            
        self.vector_store = None
        self.current_paper_id = None
        self.memory = {}
        self.pubmed_service = PubMedService()
        logger.info("ScholarRAG initialization complete")

    async def fetch_paper_data(self, paper_id: str) -> Dict:
        logger.info(f"Fetching paper data for ID: {paper_id}")
        try:
            # Use PubMedService statically, not from instance
            paper_data = await PubMedService.fetch_paper_data(paper_id)
            logger.info(f"Successfully fetched paper: {paper_data.get('title', '')[:50]}...")
            return paper_data
        except Exception as e:
            logger.error(f"Error fetching paper data: {e}")
            traceback.print_exc()
            raise

    def _create_vector_store(self, paper_data: Dict):
        logger.info("Creating vector store")
        if not paper_data:
            logger.error("Cannot create vector store: paper_data is empty")
            return
            
        logger.debug(f"Paper data: title={paper_data.get('title', '')[:30]}..., " +
                   f"authors={len(paper_data.get('authors', []))}, " +
                   f"abstract length={len(paper_data.get('abstract', ''))}")
            
        text = f"""Title: {paper_data.get('title', '')}
        Abstract: {paper_data.get('abstract', '')}
        Authors: {', '.join(paper_data.get('authors', []))}
        Year: {paper_data.get('year', '')}
        Venue: {paper_data.get('venue', '')}
        """
        
        try:
            # Use custom embeddings wrapper to avoid compatibility issues
            logger.info("Initializing embeddings model")
            embeddings_wrapper = CustomVertexEmbeddings(
                model_name="textembedding-gecko@001"
            )
            
            logger.info("Creating FAISS vector store")
            self.vector_store = FAISS.from_texts(
                texts=[text],
                embedding=embeddings_wrapper
            )
            logger.info("Vector store created successfully")
        except Exception as e:
            logger.error(f"Error creating vector store: {e}")
            traceback.print_exc()
            # Fallback to direct completion without embedding
            self.vector_store = None

    async def ask_question(self, paper_id: str, question: str) -> Dict[str, Any]:
        logger.info(f"Processing question for paper ID: {paper_id}")
        logger.info(f"Question: {question}")
        
        try:
            if question in self.memory:
                logger.info("Returning cached answer")
                return {"answer": self.memory[question], "cached": True}

            logger.info("Fetching paper data")
            paper_data = await self.fetch_paper_data(paper_id)
            
            if paper_id != self.current_paper_id:
                logger.info(f"New paper ID detected. Old: {self.current_paper_id}, New: {paper_id}")
                logger.info("Creating vector store")
                self._create_vector_store(paper_data)
                self.current_paper_id = paper_id

            if not self.vector_store:
                # Log warning but continue with direct query
                logger.warning("Vector store not initialized, falling back to direct query")

            logger.info("Generating prompt")
            prompt = f"""You are a medical research assistant. Based on the following paper information,
            provide detailed and accurate answers. If the information is insufficient to answer the question, 
            please clearly state so.

            Paper Information: {paper_data.get('title', '')}
            {paper_data.get('abstract', '')}

            Question: {question}"""

            logger.info("Sending request to Gemini")
            response = self.model.generate_content(prompt)
            answer = response.text
            logger.info(f"Received response of length: {len(answer)}")

            self.memory[question] = answer
            logger.info("Answer stored in memory")

            return {
                "answer": answer,
                "paper_data": paper_data,
                "cached": False
            }

        except Exception as e:
            logger.error(f"Error in ask_question: {str(e)}")
            traceback.print_exc()
            return {"error": f"Sorry, something went wrong. Error: {str(e)}"}
