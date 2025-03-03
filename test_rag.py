import asyncio
from medical_search import search_papers, ask_question
import logging

logging.basicConfig(level=logging.INFO)

async def test_rag_system():
    # First, search for a specific paper about COVID-19
    papers = await search_papers("COVID-19 vaccination effectiveness", limit=1)
    if not papers:
        print("No papers found!")
        return
    
    paper_id = papers[0]['paperId']
    print(f"\nTesting with paper: {papers[0]['title']}")
    print(f"Paper ID: {paper_id}")
    print("-" * 80)

    # Test cases with different types of questions
    test_questions = [
        # Basic information extraction
        "What is the main topic of this paper?",
        "What are the key findings of this study?",
        
        # Specific detail questions
        "What methodology was used in this study?",
        "What were the limitations mentioned in the study?",
        
        # Analysis questions
        "How does this study compare to previous research?",
        "What are the implications of these findings for clinical practice?",
        
        # Edge cases
        "Can you summarize the statistical analysis used?",
        "What future research directions are suggested?",
        
        # Questions about specific sections
        "What does the abstract say about the results?",
        "What were the inclusion criteria for the study?"
    ]

    for i, question in enumerate(test_questions, 1):
        print(f"\nTest Case {i}: {question}")
        print("-" * 40)
        
        try:
            response = await ask_question(paper_id, question)
            if "error" in response:
                print(f"Error: {response['error']}")
            else:
                print(f"Answer: {response['answer']}")
                if response.get("cached"):
                    print("(Retrieved from cache)")
        except Exception as e:
            print(f"Error occurred: {str(e)}")
        
        print("-" * 40)
        # Add a small delay between questions
        await asyncio.sleep(1)

if __name__ == "__main__":
    print("Starting RAG System Test...")
    print("=" * 80)
    
    try:
        asyncio.run(test_rag_system())
    except KeyboardInterrupt:
        print("\nTest interrupted by user.")
    except Exception as e:
        print(f"\nTest failed with error: {str(e)}")
    finally:
        print("\nTest completed.")
