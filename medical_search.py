from typing import Dict, Any
from langchain_core.prompts import PromptTemplate
from langchain_openai import OpenAIEmbeddings, ChatOpenAI
from langchain_community.vectorstores import Chroma
import logging
import asyncio
import aiohttp
import os
from dotenv import load_dotenv
from Bio import Entrez
import xml.etree.ElementTree as ET
import json

# Initialize logging
logging.basicConfig(level=logging.INFO)

# Load environment variables from .env file
load_dotenv()

# Get OpenAI API key
OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
if not OPENAI_API_KEY:
    logging.error("OPENAI_API_KEY not found in environment variables")
    raise ValueError("OPENAI_API_KEY not found in environment variables")

# Add PubMed configuration
Entrez.email = os.getenv("ENTREZ_EMAIL", "your-email@example.com")
Entrez.api_key = os.getenv("ENTREZ_API_KEY")  # Optional but recommended

class PubMedAPI:
    @staticmethod
    async def search_papers(query: str, limit: int = 10) -> list:
        """Search papers using PubMed E-utilities"""
        try:
            # Use asyncio to run sync Entrez functions
            loop = asyncio.get_event_loop()
            
            # Search for paper IDs
            search_results = await loop.run_in_executor(
                None,
                lambda: Entrez.read(Entrez.esearch(
                    db="pubmed",
                    term=query,
                    retmax=limit,
                    sort="relevance"
                ))
            )
            
            if not search_results["IdList"]:
                return []
            
            # Fetch paper details
            papers_xml = await loop.run_in_executor(
                None,
                lambda: Entrez.efetch(
                    db="pubmed",
                    id=search_results["IdList"],
                    rettype="xml"
                ).read()
            )
            
            # Parse XML response
            root = ET.fromstring(papers_xml)
            papers = []
            
            for article in root.findall(".//PubmedArticle"):
                try:
                    # Extract paper details
                    pmid = article.find(".//PMID").text
                    title = article.find(".//ArticleTitle").text or "No title"
                    abstract = article.find(".//Abstract/AbstractText")
                    abstract = abstract.text if abstract is not None else ""
                    
                    # Extract authors
                    authors = []
                    author_list = article.findall(".//Author")
                    for author in author_list:
                        lastname = author.find("LastName")
                        firstname = author.find("ForeName")
                        if lastname is not None and firstname is not None:
                            authors.append(f"{firstname.text} {lastname.text}")
                    
                    # Extract year
                    year_elem = article.find(".//PubDate/Year")
                    year = year_elem.text if year_elem is not None else "N/A"
                    
                    # Extract journal info
                    journal = article.find(".//Journal/Title")
                    venue = journal.text if journal is not None else "N/A"
                    
                    papers.append({
                        "paperId": pmid,
                        "title": title,
                        "abstract": abstract,
                        "authors": authors,
                        "year": year,
                        "venue": venue
                    })
                    
                except Exception as e:
                    logging.error(f"Error parsing paper: {str(e)}")
                    continue
            
            return papers
            
        except Exception as e:
            logging.error(f"Error searching PubMed: {str(e)}")
            return []

    @staticmethod
    async def fetch_paper_data(paper_id: str) -> Dict:
        """Fetch paper data from PubMed using E-utilities"""
        try:
            loop = asyncio.get_event_loop()
            
            # Fetch paper details
            paper_xml = await loop.run_in_executor(
                None,
                lambda: Entrez.efetch(
                    db="pubmed",
                    id=paper_id,
                    rettype="xml"
                ).read()
            )
            
            # Parse XML
            root = ET.fromstring(paper_xml)
            article = root.find(".//PubmedArticle")
            
            if article is None:
                raise Exception(f"Paper {paper_id} not found")
            
            # Extract paper details
            title = article.find(".//ArticleTitle").text or "No title"
            abstract = article.find(".//Abstract/AbstractText")
            abstract = abstract.text if abstract is not None else ""
            
            # Extract authors
            authors = []
            author_list = article.findall(".//Author")
            for author in author_list:
                lastname = author.find("LastName")
                firstname = author.find("ForeName")
                if lastname is not None and firstname is not None:
                    authors.append(f"{firstname.text} {lastname.text}")
            
            # Extract year and venue
            year_elem = article.find(".//PubDate/Year")
            year = year_elem.text if year_elem is not None else "N/A"
            journal = article.find(".//Journal/Title")
            venue = journal.text if journal is not None else "N/A"
            
            return {
                'title': title,
                'abstract': abstract,
                'authors': authors,
                'year': year,
                'venue': venue
            }
            
        except Exception as e:
            logging.error(f"Error fetching paper data: {str(e)}")
            raise

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
        self.pubmed_api = PubMedAPI()
        
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
        """Fetch paper data using PubMed API"""
        return await self.pubmed_api.fetch_paper_data(paper_id)

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

# Replace the existing search_papers function
async def search_papers(query: str, limit: int = 10) -> list:
    """Search papers using PubMed API"""
    return await PubMedAPI.search_papers(query, limit)

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
