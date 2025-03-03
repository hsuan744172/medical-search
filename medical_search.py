import google.generativeai as genai
from google.cloud import aiplatform
from vertexai.language_models import TextEmbeddingModel
from typing import Dict, Any
from langchain_core.prompts import PromptTemplate
import logging
import asyncio
import os
from dotenv import load_dotenv
import os
os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "/path/to/service-account.json"
# No import needed since PubMedService is defined in this file

# Initialize logging
logging.basicConfig(level=logging.INFO)

# Load environment variables from .env file
load_dotenv()

# Get Google API key from environment variables
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")
if not GOOGLE_API_KEY:
    logging.error("GOOGLE_API_KEY not found in environment variables")
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
        """Fetch paper data using PubMedService"""
        return await self.pubmed_service.fetch_paper_data(paper_id)

    def _get_embeddings(self, texts):
        """Get embeddings using Vertex AI"""
        embeddings = self.embedding_model.get_embeddings(texts)
        return [embedding.values for embedding in embeddings]

    def _create_vector_store(self, paper_data: Dict):
        """Create vector store from paper data"""
        text = f"""Title: {paper_data['title']}
        Abstract: {paper_data['abstract']}
        Authors: {', '.join(paper_data['authors'])}
        Year: {paper_data['year']}
        Venue: {paper_data['venue']}
        """
        embeddings = self._get_embeddings([text])
        # Note: You might need to modify FAISS initialization to work with Vertex AI embeddings
        self.vector_store = FAISS.from_embeddings(
            text_embeddings=list(zip([text], embeddings)),
            embedding=self.embedding_model
        )

    async def ask_question(self, paper_id: str, question: str) -> Dict[str, Any]:
        try:
            paper_data = None
            
            if question in self.memory:
                return {"answer": self.memory[question], "cached": True}

            # Fetch paper data using PubMedService
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

async def search_papers(query: str, limit: int = 10) -> list:
    """Search for papers using PubMed"""
    return await PubMedService.search_papers(query, limit)

# Initialize system
rag_system = ScholarRAG()

# Main interface
async def main_query(paper_id: str, question: str) -> Dict[str, Any]:
    response = await rag_system.ask_question(paper_id, question)
    if "error" in response:
        logging.error(response["error"])
    elif response.get("cached"):
        logging.info("Answer retrieved from cache.")
    return response

# Main interface
async def ask_question(paper_id: str, question: str) -> Dict[str, Any]:
    """Main question answering interface for compatibility"""
    return await main_query(paper_id, question)

# 確保這些是模組的公開介面
__all__ = ['ask_question', 'search_papers']
from Bio import Entrez
from typing import Dict, List
import logging
import os
from dotenv import load_dotenv

load_dotenv()
Entrez.email = os.getenv("ENTREZ_EMAIL")
if api_key := os.getenv("ENTREZ_API_KEY"):
    Entrez.api_key = api_key

class PubMedService:
    @staticmethod
    async def search_papers(query: str, max_results: int = 10) -> List[Dict]:
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
            records = Entrez.read(handle)
            handle.close()

            if not records['IdList']:
                return []

            handle = Entrez.efetch(db="pubmed", id=records['IdList'], rettype="medline", retmode="xml")
            papers = Entrez.parse(handle)
            
            results = []
            for paper in papers:
                paper_dict = {
                    'paperId': paper.get('PMID', [''])[0],
                    'title': paper.get('TI', [''])[0],
                    'abstract': paper.get('AB', [''])[0] if paper.get('AB') else '',
                    'authors': [author['LastName'] + ' ' + author['ForeName'] 
                              for author in paper.get('AU', [])],
                    'year': paper.get('DP', [''])[0][:4],
                    'venue': paper.get('JT', [''])[0] if paper.get('JT') else paper.get('TA', [''])[0],
                }
                results.append(paper_dict)
            
            handle.close()
            return results
        except Exception as e:
            logging.error(f"Error searching PubMed: {str(e)}")
            return []

    @staticmethod
    async def fetch_paper_data(paper_id: str) -> Dict:
        try:
            handle = Entrez.efetch(db="pubmed", id=paper_id, rettype="medline", retmode="xml")
            papers = list(Entrez.parse(handle))
            if not papers:
                raise Exception(f"No paper found with ID {paper_id}")
            
            paper = papers[0]
            paper_data = {
                'title': paper.get('TI', [''])[0],
                'abstract': paper.get('AB', [''])[0] if paper.get('AB') else '',
                'authors': [author['LastName'] + ' ' + author['ForeName'] 
                          for author in paper.get('AU', [])],
                'year': paper.get('DP', [''])[0][:4],
                'venue': paper.get('JT', [''])[0] if paper.get('JT') else paper.get('TA', [''])[0],
                'keywords': paper.get('MH', [])
            }
            handle.close()
            return paper_data
            
        except Exception as e:
            logging.error(f"Error fetching paper data: {str(e)}")
            raise
