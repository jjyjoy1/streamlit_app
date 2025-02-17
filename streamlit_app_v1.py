import streamlit as st

create_page = st.Page("second_step_process_streamlit_prod_v3.py", title="Soham_Tool", icon=":material/add_circle:")
create_page2 = st.Page("first_step_process_streamlit_pord_v2.py", title="BD_Tools", icon=":material/add_circle:")

pg = st.navigation([create_page, create_page2])
st.set_page_config(page_title="Centralized Data manager", page_icon=":material/edit:")
pg.run()


