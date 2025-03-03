# Medical Paper Q&A Assistant

A bilingual (English/Chinese) medical research assistant that helps you search and understand medical papers using Google's Gemini AI, PubMed API, and Streamlit.

## Features

- 🔍 **Paper Search**: Search medical papers through PubMed's extensive database
- 💬 **Interactive Q&A**: Ask questions about specific papers and get AI-powered responses
- 🌐 **Bilingual Interface**: Switch between English and Chinese interfaces
- 📚 **Paper Details**: View paper abstracts, authors, publication year, and venue
- 💾 **Memory Cache**: Caches responses for repeated questions to improve performance

## Requirements

- Python 3.11 or higher
- Google Cloud Project with Gemini API enabled
- PubMed API access (email required, API key recommended)

## Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/medical-search.git
   cd medical-search
   ```

2. Install dependencies using Poetry:
   ```bash
   poetry install
   ```

3. Set up your environment variables in a `.env` file:
   ```
   GOOGLE_API_KEY=your_google_api_key
   ENTREZ_EMAIL=your_email@example.com
   ENTREZ_API_KEY=your_pubmed_api_key  # Optional but recommended
   ```

4. Set up Google Cloud credentials:
   - Download your service account key JSON file
   - Update the path in `medical_search.py`:
     ```python
     os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "/path/to/service-account.json"
     ```

5. Run the application:
   ```bash
   poetry run streamlit run app.py
   ```

## Usage

1. The interface will open in your default web browser
2. Use the sidebar to search for medical papers using keywords
3. Click on a paper from the search results to select it
4. Ask questions about the selected paper in the chat interface
5. Toggle between English and Chinese using the language button (🌐) in the header

## Technical Stack

- **Frontend**: Streamlit
- **AI Model**: Google Gemini Pro
- **Embeddings**: Vertex AI TextEmbedding
- **Data Source**: PubMed API via Biopython
- **Vector Store**: FAISS
- **Dependencies Management**: Poetry

## Notes

- The system uses PubMed's API, which requires an email address for requests
- Using a PubMed API key is recommended for higher rate limits
- The chat interface retains conversation history during the session
- Paper details are cached to improve performance

## License

MIT License

