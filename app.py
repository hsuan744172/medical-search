import streamlit as st
import asyncio
from medical_search import main_query as ask_question, search_papers

# 定義翻譯文本
TRANSLATIONS = {
    'en': {
        'title': 'Medical Paper Q&A Assistant',
        'search_header': 'Search Papers',
        'search_placeholder': 'e.g. COVID-19 treatment',
        'search_button': 'Search',
        'searching': 'Searching papers...',
        'search_results': 'Search Results',
        'selected_paper': 'Selected Paper',
        'show_abstract': 'Show Abstract',
        'paper_details': 'Paper Details',
        'ask_placeholder': 'Ask a question about the paper',
        'thinking': 'Thinking...',
        'search_hint': '👈 Search for papers using keywords in the sidebar',
        'title_label': 'Title:',
        'year_label': 'Year:',
        'authors_label': 'Authors:',
        'venue_label': 'Venue:'
    },
    'zh': {
        'title': '醫學論文問答助手',
        'search_header': '搜尋論文',
        'search_placeholder': '例如: COVID-19 治療',
        'search_button': '搜尋',
        'searching': '搜尋論文中...',
        'search_results': '搜尋結果',
        'selected_paper': '已選擇的論文',
        'show_abstract': '顯示摘要',
        'paper_details': '論文詳情',
        'ask_placeholder': '請輸入關於論文的問題',
        'thinking': '思考中...',
        'search_hint': '👈 在側邊欄使用關鍵字搜尋論文',
        'title_label': '標題:',
        'year_label': '年份:',
        'authors_label': '作者:',
        'venue_label': '發表於:'
    }
}

def init_session_state():
    """Initialize session state variables"""
    if 'paper_id' not in st.session_state:
        st.session_state.paper_id = ''
    if 'chat_history' not in st.session_state:
        st.session_state.chat_history = []
    if 'search_results' not in st.session_state:
        st.session_state.search_results = []
    if 'selected_paper' not in st.session_state:
        st.session_state.selected_paper = None
    if 'language' not in st.session_state:
        st.session_state.language = 'en'

def get_text(key):
    """Get translated text based on current language"""
    return TRANSLATIONS[st.session_state.language][key]

def main():
    init_session_state()
    
    # Language selector in header
    col1, col2 = st.columns([6, 1])
    with col1:
        st.title(get_text('title'))
    with col2:
        if st.button('🌐 ' + ('EN' if st.session_state.language == 'zh' else '中')):
            st.session_state.language = 'en' if st.session_state.language == 'zh' else 'zh'
            st.rerun()

    # Search interface
    with st.sidebar:
        st.header(get_text('search_header'))
        search_query = st.text_input(
            get_text('search_header') + ":",
            placeholder=get_text('search_placeholder')
        )
        
        if st.button(get_text('search_button')):
            with st.spinner(get_text('searching')):
                st.session_state.search_results = asyncio.run(search_papers(search_query))
        
        # Display search results
        if st.session_state.search_results:
            st.subheader(get_text('search_results'))
            for paper in st.session_state.search_results:
                if st.button(
                    f"{paper['title'][:100]}...",
                    key=paper['paperId'],
                    help=f"{get_text('authors_label')} {', '.join(paper['authors'])}\n{get_text('year_label')} {paper['year']}"
                ):
                    st.session_state.paper_id = paper['paperId']
                    st.session_state.selected_paper = paper
                    st.session_state.chat_history = []
                    st.rerun()

        # Show selected paper details
        if st.session_state.selected_paper:
            st.divider()
            st.subheader(get_text('selected_paper'))
            paper = st.session_state.selected_paper
            st.markdown(f"**{get_text('title_label')}** {paper['title']}")
            st.markdown(f"**{get_text('year_label')}** {paper['year']}")
            st.markdown(f"**{get_text('authors_label')}** {', '.join(paper['authors'])}")
            if paper.get('abstract'):
                with st.expander(get_text('show_abstract')):
                    st.write(paper['abstract'])

    # Chat interface
    if st.session_state.paper_id:
        # Display chat history
        for i, (question, answer) in enumerate(st.session_state.chat_history):
            with st.chat_message("user"):
                st.write(question)
            with st.chat_message("assistant"):
                st.write(answer)

        # Question input
        if question := st.chat_input(get_text('ask_placeholder')):
            with st.chat_message("user"):
                st.write(question)

            with st.chat_message("assistant"):
                with st.spinner(get_text('thinking')):
                    response = asyncio.run(ask_question(st.session_state.paper_id, question))
                    
                    if "error" in response:
                        st.error(response["error"])
                    else:
                        answer = response["answer"]
                        st.write(answer)
                        st.session_state.chat_history.append((question, answer))

                        if len(st.session_state.chat_history) == 1 and "paper_data" in response:
                            with st.sidebar:
                                st.subheader(get_text('paper_details'))
                                paper_data = response["paper_data"]
                                st.write(f"**{get_text('title_label')}** {paper_data['title']}")
                                st.write(f"**{get_text('year_label')}** {paper_data['year']}")
                                st.write(f"**{get_text('authors_label')}** {', '.join(paper_data['authors'])}")
                                st.write(f"**{get_text('venue_label')}** {paper_data['venue']}")
    else:
        st.info(get_text('search_hint'))

if __name__ == "__main__":
    main()
