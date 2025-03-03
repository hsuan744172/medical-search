from typing import Dict, Any
from langchain_core.prompts import PromptTemplate
from langchain_openai import OpenAIEmbeddings, ChatOpenAI
from langchain_community.vectorstores import Chroma
import logging
import asyncio
import aiohttp
import os
from dotenv import load_dotenv

# Initialize logging
logging.basicConfig(level=logging.INFO)

# Load environment variables from .env file
load_dotenv()

# Get OpenAI API key
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
if not OPENAI_API_KEY:
    logging.error("OPENAI_API_KEY not found in environment variables")
    raise ValueError("OPENAI_API_KEY not found in environment variables")

class ScholarRAG:
    def __init__(self):
        self.embeddings = self._init_embeddings()
        self.llm = self._init_llm()
        self.vector_store = None
        self.current_paper_id = None
        self.memory = {}
        self.api_url = "https://api.semanticscholar.org/v1/paper"
        self.search_api_url = "https://api.semanticscholar.org/graph/v1/paper/search"
        self.persist_directory = "chroma_db"  # Add persist directory for Chroma
        
    def _init_embeddings(self):
        """Initialize OpenAI embedding model"""
        return OpenAIEmbeddings(
            model="text-embedding-3-small"
        )
    
    def _init_llm(self):
        """Initialize OpenAI GPT-4 Turbo"""
        return ChatOpenAI(
            model="gpt-4-0125-preview",
            temperature=0.7
        )

    async def fetch_paper_data(self, paper_id: str) -> Dict:
        """Fetch paper data asynchronously"""
        headers = {"User-Agent": "Mozilla/5.0"}
        url = f"{self.api_url}/{paper_id}"
        async with aiohttp.ClientSession() as session:
            try:
                async with session.get(url, headers=headers) as response:
                    if response.status != 200:
                        raise Exception(f"Failed to fetch paper {paper_id}, status code {response.status}")
                    data = await response.json()
                    return {
                        'title': data.get('title', ''),
                        'abstract': data.get('abstract', ''),
                        'authors': [author.get('name') for author in data.get('authors', [])],
                        'year': data.get('year'),
                        'venue': data.get('venue'),
                        'citations': data.get('citations', [])
                    }
            except Exception as e:
                logging.error(f"Error fetching paper data: {str(e)}")
                raise

    def _create_vector_store(self, paper_data: Dict):
        """Create vector store from paper data"""
        text = f"""Title: {paper_data['title']}
        Abstract: {paper_data['abstract']}
        Authors: {', '.join(paper_data['authors'])}
        Year: {paper_data['year']}
        Venue: {paper_data['venue']}
        """
        self.vector_store = Chroma.from_texts(
            texts=[text],
            embedding=self.embeddings,
            persist_directory=self.persist_directory,
            collection_name=f"paper_{self.current_paper_id}"
        )

    async def ask_question(self, paper_id: str, question: str) -> Dict[str, Any]:
        try:
            paper_data = None  # Initialize paper_data outside try block
            
            if question in self.memory:
                return {"answer": self.memory[question], "cached": True}

            if paper_id != self.current_paper_id:
                paper_data = await self.fetch_paper_data(paper_id)
                self._create_vector_store(paper_data)
                self.current_paper_id = paper_id
            else:
                # Fetch paper data even if it's the same paper
                paper_data = await self.fetch_paper_data(paper_id)

            if not self.vector_store:
                raise Exception("Vector store not initialized")

            prompt_template = """You are a medical research assistant. Based on the following paper information,
            provide detailed and accurate answers. If the information is insufficient to answer the question, please clearly state so.

            Paper Information: {context}
            Question: {question}"""
            
            PROMPT = PromptTemplate(
                template=prompt_template,
                input_variables=["context", "question"]
            )

            docs = self.vector_store.similarity_search(question, k=1)
            context = "\n".join(doc.page_content for doc in docs)
            formatted_prompt = PROMPT.format(context=context, question=question)
            
            answer = self.llm.invoke(formatted_prompt).content

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
    """Search for papers using Semantic Scholar API"""
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"
    }
    
    params = {
        "query": query,
        "limit": limit,
        "fields": "title,abstract,authors,year,venue,paperId"
    }
    
    try:
        async with aiohttp.ClientSession() as session:
            async with session.get(
                "https://api.semanticscholar.org/graph/v1/paper/search",
                params=params,
                headers=headers
            ) as response:
                if response.status != 200:
                    raise Exception("Failed to search papers")
                    
                data = await response.json()
                papers = data.get("data", [])
                
                # Format the results
                return [{
                    "title": paper.get("title", ""),
                    "abstract": paper.get("abstract", ""),
                    "authors": [author.get("name", "") for author in paper.get("authors", [])],
                    "year": paper.get("year"),
                    "venue": paper.get("venue", ""),
                    "paperId": paper.get("paperId", "")
                } for paper in papers]
                
    except Exception as e:
        print(f"Error searching papers: {str(e)}")
        return []

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
