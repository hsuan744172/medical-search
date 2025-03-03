from typing import Dict, Any
import logging
from .models.scholar_rag import ScholarRAG
from .services.pubmed_service import PubMedService

# Initialize system
rag_system = ScholarRAG()

async def search_papers(query: str, limit: int = 10) -> list:
    """Search for papers using PubMed"""
    return await PubMedService.search_papers(query, limit)

async def ask_question(paper_id: str, question: str) -> Dict[str, Any]:
    """Question answering interface"""
    return await rag_system.ask_question(paper_id, question)

__all__ = ['ask_question', 'search_papers']
