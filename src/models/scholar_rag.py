import google.generativeai as genai
from google.cloud import aiplatform
from vertexai.language_models import TextEmbeddingModel
from langchain_community.vectorstores import FAISS
from langchain_community.embeddings import VertexAIEmbeddings
from typing import Dict, Any
import logging
import os
from dotenv import load_dotenv
from ..services.pubmed_service import PubMedService

# Initialize logging
logging.basicConfig(level=logging.INFO)

# Load environment variables
load_dotenv()

# Configure Google credentials
os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "/path/to/service-account.json"
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")
if not GOOGLE_API_KEY:
    raise ValueError("GOOGLE_API_KEY not found in environment variables")

# Configure Gemini
genai.configure(api_key=GOOGLE_API_KEY)

class ScholarRAG:
    def __init__(self):
        self.model = genai.GenerativeModel('gemini-pro')
        self.embedding_model = TextEmbeddingModel.from_pretrained("textembedding-gecko@001")
        self.vector_store = None
        self.current_paper_id = None
        self.memory = {}
        self.pubmed_service = PubMedService()

    async def fetch_paper_data(self, paper_id: str) -> Dict:
        return await self.pubmed_service.fetch_paper_data(paper_id)

    def _get_embeddings(self, texts):
        embeddings = self.embedding_model.get_embeddings(texts)
        return [embedding.values for embedding in embeddings]

    def _create_vector_store(self, paper_data: Dict):
        text = f"""Title: {paper_data['title']}
        Abstract: {paper_data['abstract']}
        Authors: {', '.join(paper_data['authors'])}
        Year: {paper_data['year']}
        Venue: {paper_data['venue']}
        """
        
        embeddings_wrapper = VertexAIEmbeddings(
            model_name="textembedding-gecko@001"
        )
        
        self.vector_store = FAISS.from_texts(
            texts=[text],
            embedding=embeddings_wrapper
        )

    async def ask_question(self, paper_id: str, question: str) -> Dict[str, Any]:
        try:
            if question in self.memory:
                return {"answer": self.memory[question], "cached": True}

            paper_data = await self.fetch_paper_data(paper_id)
            
            if paper_id != self.current_paper_id:
                self._create_vector_store(paper_data)
                self.current_paper_id = paper_id

            if not self.vector_store:
                raise Exception("Vector store not initialized")

            prompt = f"""You are a medical research assistant. Based on the following paper information,
            provide detailed and accurate answers. If the information is insufficient to answer the question, 
            please clearly state so.

            Paper Information: {paper_data['title']}
            {paper_data['abstract']}

            Question: {question}"""

            response = self.model.generate_content(prompt)
            answer = response.text

            self.memory[question] = answer

            return {
                "answer": answer,
                "paper_data": paper_data,
                "cached": False
            }

        except Exception as e:
            logging.error(f"Error in ask_question: {str(e)}")
            return {"error": f"Sorry, something went wrong. Error: {str(e)}"}
