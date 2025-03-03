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
