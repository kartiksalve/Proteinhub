import streamlit as st
import requests
import networkx as nx
import plotly.graph_objects as go
# import base64 # No longer needed for background
# import os # Not used
import io

# -------------- Background Styling (Removed background image, kept container styling) --------------
def style_app_container():
    css = f"""
    <style>
    .block-container {{
        background-color: rgba(0, 0, 0, 0.6); /* Kept for content readability if desired */
        /* You can adjust or remove this background-color if you want a plain white background for the container */
        /* background-color: white; */ /* Example for a white container */
        padding: 2rem 3rem;
        border-radius: 1rem;
        /* backdrop-filter: blur(8px); */ /* Removed as it was tied to the background image effect */
        /* -webkit-backdrop-filter: blur(8px); */ /* Removed */
        color: #333; /* Changed to a more standard text color for better readability without dark background */
        /* color: #f0f0f0; */ /* Original color for dark background */
        max-width: 1000px;
        margin: auto;
    }}
    /* Optional: Style for the main app body if you want a specific color */
    .stApp {{
        /* background-color: #f0f2f6; */ /* Example: Light grey background for the whole app */
    }}
    </style>
    """
    st.markdown(css, unsafe_allow_html=True)

style_app_container() # Apply container styling

# -------------- Constants --------------
STRING_API_URL = "https://string-db.org/api"
STRING_OUTPUT_FORMAT = "json"
STRING_METHOD = "network"

# -------------- Functions --------------
def get_string_interactions(uniprot_id, species=9606, min_score=0.4):
    params = {
        "identifiers": uniprot_id,
        "species": species,
        "caller_identity": "streamlit_app",
        "required_score": int(min_score * 1000)
    }
    url = f"{STRING_API_URL}/{STRING_OUTPUT_FORMAT}/{STRING_METHOD}"
    try:
        response = requests.post(url, data=params)
        response.raise_for_status()  # Raise an exception for HTTP errors
        return response.json()
    except requests.exceptions.RequestException as e:
        st.error(f"Error fetching data from STRING DB: {e}")
        return None
    except ValueError: # Catches JSON decoding errors
        st.error("Error decoding JSON response from STRING DB. The API might have returned an unexpected format.")
        if response:
            st.error(f"Response content: {response.text[:500]}...") # Show part of the response
        return None


def build_network(data):
    G = nx.DiGraph() # Using DiGraph as in original, though STRING interactions are often undirected.
                    # If undirected is preferred, use nx.Graph()
    if not data or not isinstance(data, list): # Add check for empty or malformed data
        return G
        
    for interaction in data:
        # Ensure keys exist to prevent KeyErrors
        p1 = interaction.get('preferredName_A')
        p2 = interaction.get('preferredName_B')
        score = interaction.get('score')
        
        if p1 and p2 and score is not None: # Ensure all required fields are present
            G.add_edge(p1, p2, weight=float(score)) # Ensure score is float
    return G

def find_hub_genes(G, top_n=5):
    if G.number_of_nodes() == 0:
        return []
    degree_dict = dict(G.degree())
    sorted_genes = sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)
    return [gene for gene, _ in sorted_genes[:top_n]]

def create_graph_figure(G, hub_genes):
    if G.number_of_nodes() == 0: # Handle empty graph case for layout
        pos = {}
    else:
        pos = nx.spring_layout(G, seed=42)
    
    degrees = dict(G.degree())

    edge_x, edge_y = [], []
    for src, dst in G.edges():
        if src in pos and dst in pos: # Ensure nodes are in pos dict
            x0, y0 = pos[src]
            x1, y1 = pos[dst]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=1.5, color='gray'),
        hoverinfo='none',
        mode='lines'
    )

    node_x, node_y, node_size, node_color, node_text_labels, node_hover_text = [], [], [], [], [], []
    for node in G.nodes():
        if node in pos: # Ensure node is in pos dict
            x, y = pos[node]
            degree = degrees.get(node, 0) # Use .get for safety
            node_x.append(x)
            node_y.append(y)
            node_size.append(15 + degree * 2)
            node_color.append('red' if node in hub_genes else 'royalblue')
            node_text_labels.append(node) # Text displayed on the node
            node_hover_text.append(f"{node}<br>Degree: {degree}") # Text for hover

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        text=node_text_labels, # Text elements to display on the graph
        textposition="middle center",
        textfont=dict(color='white', size=10), # Text on nodes (inside markers if small enough)
        marker=dict(size=node_size, color=node_color, line=dict(width=1, color='white')),
        hoverinfo='text', # Use 'text' and supply hover_text
        hovertext=node_hover_text # Custom hover text
    )

    # The main fix is here: titlefont -> title_font
    layout = go.Layout(
        title="Protein Interaction Network",
        title_font=dict(color='#333'), # Changed from titlefont
        showlegend=False,
        hovermode='closest',
        margin=dict(b=20, l=20, r=20, t=40),
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        plot_bgcolor='rgba(255,255,255,1)', # White plot background
        paper_bgcolor='rgba(255,255,255,1)', # White paper background
    )

    return go.Figure(data=[edge_trace, node_trace], layout=layout)

# -------------- Streamlit UI --------------
st.set_page_config(page_title="Prot'n'Hub", layout="wide") # Optional: use wide layout

st.title("🧬 Prot'n'Hub – Protein Interaction & Hub Gene Explorer")

tabs = st.tabs(["Home", "About"])
with tabs[0]:
    st.header("Explore Protein Network")

    user_input = st.text_area("Enter Protein Name or UniProt ID(s)", height=100, 
                              help="Enter one or more UniProt IDs or protein names, separated by newlines or commas.")

    st.subheader("Species Selection")
    species_dict = {
        "Human (Homo sapiens)": 9606,
        "Mouse (Mus musculus)": 10090,
        "Rat (Rattus norvegicus)": 10116,
        "Zebrafish (Danio rerio)": 7955,
        "Fruit fly (Drosophila melanogaster)": 7227,
        "Yeast (Saccharomyces cerevisiae)": 4932, # Added yeast
        "E. coli (Escherichia coli K12)": 83333, # Added E. coli
        "Custom (enter manually)": None,
    }
    selected_species_name = st.selectbox("Choose species", list(species_dict.keys()), index=0)
    if selected_species_name == "Custom (enter manually)":
        species = st.number_input("Enter NCBI Taxonomy ID", value=9606, min_value=0)
    else:
        species = species_dict[selected_species_name]

    score_threshold = st.slider("Interaction Score Threshold", 0.0, 1.0, 0.4, 0.05,
                                help="Minimum confidence score for interactions (0.0-1.0). Higher values mean more stringent, fewer interactions.")

    if st.button("🔍 Analyze", type="primary"):
        with st.spinner("Fetching and analyzing data..."):
            # Process multiple inputs (comma or newline separated)
            raw_input = user_input.strip()
            if not raw_input:
                st.warning("Please enter at least one Protein Name or UniProt ID.")
                st.stop()
            
            # Split by comma or newline, strip whitespace, filter empty strings
            identifiers_list = [ident.strip() for ident_group in raw_input.split(',') 
                                for ident in ident_group.splitlines() if ident.strip()]
            
            if not identifiers_list:
                st.warning("No valid identifiers found after processing input.")
                st.stop()

            # Join identifiers with '%0d' for STRING API (newline character)
            uniprot_ids_for_api = "%0d".join(identifiers_list)

            string_data = get_string_interactions(uniprot_ids_for_api, species, score_threshold)
            
            if string_data is None: # get_string_interactions now returns None on error and logs it
                # Error message already shown by get_string_interactions
                st.stop()
            if not string_data: # Empty list from API (valid response, but no interactions)
                 st.info(f"ℹ️ No interaction data found for '{', '.join(identifiers_list)}' with the current settings. "
                         "This could be due to the protein(s) not having known interactions in STRING, "
                         "the species, or the score threshold. Try a lower score or check your input.")
                 st.stop()


            G = build_network(string_data)
            if G.number_of_nodes() == 0:
                st.warning("ℹ️ No network could be built with the given parameters. "
                           "This might happen if the provided identifiers were not found or had no interactions "
                           "above the selected score threshold. Try a lower score threshold or check your identifiers.")
                st.stop()

            hub_genes = find_hub_genes(G)
            if hub_genes:
                st.success(f"Top Hub Genes: {', '.join(hub_genes)}")
            else:
                st.info("No distinct hub genes found based on the current network (e.g., all nodes have similar degrees or network is too small).")


            fig = create_graph_figure(G, hub_genes)
            st.plotly_chart(fig, use_container_width=True)

            # Save to PNG and offer download
            try:
                buf = io.BytesIO()
                fig.write_image(buf, format="png", width=1000, height=800, engine="kaleido") # Ensure kaleido is installed
                st.download_button(
                    label="📥 Download Network as PNG",
                    data=buf.getvalue(),
                    file_name=f"{identifiers_list[0]}_network.png", # Use first identifier for filename
                    mime="image/png"
                )
            except Exception as e:
                st.error(f"Could not generate PNG for download. Make sure 'kaleido' is installed (`pip install kaleido`). Error: {e}")


            with st.expander("📊 Network Analysis Results", expanded=True):
                st.write(f"⭐ **Nodes (Proteins)**: {G.number_of_nodes()}")
                # st.write(f"List of Nodes: {', '.join(list(G.nodes()))}") # Can be too long
                st.write(f"⭐ **Edges (Interactions)**: {G.number_of_edges()}")
                
                degree_dict = dict(G.degree())
                if degree_dict:
                    st.write("⭐ **Node Degrees (Top 10 or all if fewer):**")
                    sorted_degrees = sorted(degree_dict.items(), key=lambda item: item[1], reverse=True)
                    for node, degree in sorted_degrees[:10]:
                        st.write(f"- {node}: {degree}")
                    if len(sorted_degrees) > 10:
                        st.write(f"... and {len(sorted_degrees)-10} more.")

                    main_hub_node, main_hub_degree = sorted_degrees[0]
                    st.info(f"🏆 **Main Hub Gene (Highest Degree):** {main_hub_node} (Degree: {main_hub_degree})")
                else:
                    st.info("No nodes to determine degrees or main hub gene.")

with tabs[1]:
    st.header("About Prot'n'Hub")
    st.markdown("""
    Prot'n'Hub is a Streamlit-based interactive application for exploring protein-protein interaction networks.

    **Features:**
    - Input UniProt ID(s) or protein name(s)
    - Visualize interaction graphs
    - Detect top hub genes
    - Interactive, styled Plotly graphs
    - Customizable species and score thresholds
    - Downloadable graph as PNG image

    **Powered By:**
    - STRING API
    - UniProt API (Implicitly via STRING for name resolution)
    - NetworkX, Plotly, and Streamlit

    ---

    ### 🧑‍🏫 Quick Guide

      Prot'n'Hub is designed to be simple and informative. Here's a quick guide:

    **🔹 Explore Protein Network:**  
    Enter one or more *protein names* (like "TP53", "EGFR") or *UniProt IDs* (like "P04637", "P00533"), separated by commas or newlines. The app uses these to search for protein-protein interactions from the STRING database.

    **🔹 Species Selection:**  
    Choose the species your protein(s) belong to. For example:
    - Human = Homo sapiens
    - Mouse = Mus musculus  
    If your species isn't listed, choose **Custom** and enter its **NCBI Taxonomy ID** (a unique number for each species).

    **🔹 Interaction Score Threshold:**  
    This sets the minimum confidence score for the interactions (0.0 to 1.0).  
    - Lower values (e.g. 0.2) = more connections, but potentially lower reliability  
    - Higher values (e.g. 0.7+) = fewer connections, but higher reliability  
    Default is **0.4**, which often provides a good balance.

    **🔹 Nodes:**  
    Each circle in the graph is a *protein*. The number of nodes shows how many proteins are in your interaction network.

    **🔹 Edges:**  
    Lines connecting the nodes. Each edge represents an interaction between two proteins.

    **🔹 Node Degrees:**  
    This shows how many connections (edges) each protein (node) has.  
    Proteins with high degree values are often *hub genes*—key proteins that interact with many others.

    **🔹 Main Hub Gene:**  
    The protein with the most connections in your network. These are often biologically important and can be potential targets for further research.

    **🔹 Download Graph as PNG:**  
    Click this button to save the visual interaction network as a PNG image for reports or presentations.

    ---

    ### 📄 Acknowledgement
      I would like to express my sincere gratitude to the following resources and individuals who contributed to the development of this application:

    **STRING Database**: For providing comprehensive protein-protein interaction data, which forms the core of this application's functionality.  
    **UniProt Knowledgebase**: For enabling the mapping of protein sequences to UniProt IDs, a crucial step in processing user-provided sequence input.  
    **NetworkX**: For the powerful tools used in network analysis and manipulation, allowing for the construction and processing of protein interaction graphs.  
    **Plotly**: For the creation of interactive and visually appealing network visualizations, enhancing the user experience.  
    **Streamlit**: For providing a user-friendly framework for building and deploying the web application.  

    I extend my deepest gratitude to **Dr. Kushagra Kashyap**, for his invaluable guidance, support, and expertise throughout this project. His insights and encouragement were instrumental  in shaping the application and overcoming the challenges encountered during its development.

    ---

    **Developer Information:**

    Hi, I'm Kartik Parag Salve, the creator of Prot'n'Hub.
    As someone currently pursuing my Master's in Bioinformatics at Deccan Education Society Pune University,
    I've always been fascinated by the power of protein-protein interaction networks.
    That's why I built Prot'n'Hub which is a Streamlit app that makes exploring these interactions intuitive and engaging.

    I'm a big fan of the rich data available through the STRING and UniProt APIs,
    and the amazing visualization capabilities of NetworkX and Plotly.
    Prot'n'Hub lets you:
    - Input UniProt IDs or protein names
    - Visualize interaction graphs
    - Find key hub genes
    - Customize species and score thresholds
    - Download the graph for reports or presentations

    It's been exciting to bring these tools together in a user-friendly application,
    and I hope Prot'n'Hub proves useful for others exploring the fascinating world of protein interactions.

    **Contact:** 3522411023@despu.edu.in
    """)
