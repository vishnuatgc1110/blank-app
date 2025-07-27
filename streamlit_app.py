import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from datetime import datetime
import subprocess
from tempfile import NamedTemporaryFile
import os
import json  # Import json for JSON export

# ───────────────────────────────────────────────────────────────
# PATH CONFIGURATION(Windows)
# ───────────────────────────────────────────────────────────────
BLAST_DB = r"/workspaces/blank-app/vir_db"
BLASTN_PATH = r"/workspaces/ncbi-blast-2.17.0+/bin/blastn"

# Page configuration
st.set_page_config(
    page_title="S. maltophilia Virulence Gene Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
        padding: 2rem;
        border-radius: 10px;
        margin-bottom: 2rem;
        color: white;
        text-align: center;
    }
    .metric-card {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 1rem;
        border-radius: 10px;
        color: white;
        text-align: center;
        margin: 0.5rem 0;
    }
    .sequence-card {
        border-left: 4px solid #3498db;
        padding: 1rem;
        background: #f8f9fa;
        border-radius: 0 8px 8px 0;
        margin: 1rem 0;
    }
    .blast-match {
        border: 1px solid #e9ecef;
        border-radius: 8px;
        padding: 1rem;
        background: #f8f9fa;
        margin: 0.5rem 0;
    }
    .virulence-factor {
        background: #28a745;
        color: white;
        padding: 0.25rem 0.5rem;
        border-radius: 15px;
        font-size: 0.8rem;
        margin: 0.2rem;
        display: inline-block;
    }
</style>
""", unsafe_allow_html=True)

class StenotrophomonasAnalyzer:
    def __init__(self):
        self.virulence_genes = {
            'smlt': 'Core virulence factors',
            'rml': 'Lipopolysaccharide biosynthesis',
            'sme': 'Beta-lactamase production', 
            'pil': 'Pilus assembly',
            'tol': 'Toluene degradation',
            'xcp': 'Type II secretion',
            'hrp': 'Type III secretion',
            'flagellar': 'Motility and chemotaxis',
            'biofilm': 'Biofilm formation',
            'iron': 'Iron acquisition',
            'adhesin': 'Cell adhesion',
            'protease': 'Protease activity',
            'hemolysis': 'Hemolytic activity',
            'resistance': 'Antibiotic resistance'
        }
         # Mapping of gene names to subject IDs
        self.subject_id_mapping = {
            'acr': '1.acr',
            'bioA': '3.bioA',
            'emr': '4.emr',
            'CPR': '5.CPR',
            'feoB': '6.feoB',
            'pilU': '7.pilU'
        }
    def parse_fasta(self, file_content):
        """Parse FASTA format sequences"""
        # Implementation remains the same as your original
        sequences = []
        lines = file_content.decode('utf-8').split('\n')
        current_sequence = None
        
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(current_sequence)
                current_sequence = {
                    'header': line[1:],
                    'sequence': '',
                    'length': 0,
                    'gc_content': 0
                }
            elif line and current_sequence:
                current_sequence['sequence'] += line.upper()
        
        if current_sequence:
            sequences.append(current_sequence)
        
        for seq in sequences:
            seq['length'] = len(seq['sequence'])
            seq['gc_content'] = self.calculate_gc_content(seq['sequence'])
            seq['virulence_factors'] = self.identify_virulence_factors(seq['header'], seq['sequence'])
        
        return sequences
    
    def calculate_gc_content(self, sequence):
        """Calculate GC content percentage"""
        if not sequence:
            return 0
        gc_count = sequence.count('G') + sequence.count('C')
        return round((gc_count / len(sequence)) * 100, 2)
    
    def identify_virulence_factors(self, header, sequence):
        """Identify potential virulence factors"""
        factors = []
        header_lower = header.lower()
        
        for gene, description in self.virulence_genes.items():
            if gene in header_lower:
                factors.append({
                    'gene': gene,
                    'description': description,
                    'confidence': 'High'
                })
        
        virulence_patterns = {
            'virulence': 'General virulence factor',
            'pathogen': 'Pathogenicity factor',
            'toxin': 'Toxin production',
            'secretion': 'Secretion system',
            'efflux': 'Efflux pump'
        }
        
        for pattern, desc in virulence_patterns.items():
            if pattern in header_lower and not any(f['description'] == desc for f in factors):
                factors.append({
                    'gene': pattern,
                    'description': desc,
                    'confidence': 'Medium'
                })
        
        return factors
    
    class StenotrophomonasAnalyzer:
        def __init__(self):
          self.virulence_genes = {
            'smlt': 'Core virulence factors',
            'rml': 'Lipopolysaccharide biosynthesis',
            'sme': 'Beta-lactamase production', 
            'pil': 'Pilus assembly',
            'tol': 'Toluene degradation',
            'xcp': 'Type II secretion',
            'hrp': 'Type III secretion',
            'flagellar': 'Motility and chemotaxis',
            'biofilm': 'Biofilm formation',
            'iron': 'Iron acquisition',
            'adhesin': 'Cell adhesion',
            'protease': 'Protease activity',
            'hemolysis': 'Hemolytic activity',
            'resistance': 'Antibiotic resistance'
        }
        
        # Mapping of gene names to subject IDs
          self.subject_id_mapping = {
              'acr': '1.acr',
              'bioA': '3.bioA',
              'emr': '4.emr',
              'CPR': '5.CPR',
              'feoB': '6.feoB',
              'pilU': '7.pilU'
        }
    
    def simulate_blast_analysis(self, sequences):
        """Simulate BLAST analysis results with improved confidence levels"""
        blast_results = []
        
        for i, seq in enumerate(sequences):
            seq_id = f"Sequence_{i+1}"
            
            # Generate matches for all virulence factors (HIGH confidence)
            for vf in seq['virulence_factors']:
                gene = vf['gene']
                subject_id = self.subject_id_mapping.get(gene, f"{gene}_reference")  # Use mapping
                
                blast_results.append({
                    'query_id': seq_id,
                    'subject_id': subject_id,
                    'identity': float(np.random.uniform(90, 99)),
                    'alignment_length': int(seq['length'] * np.random.uniform(0.8, 0.95)),
                    'e_value': float(10 ** (-np.random.uniform(100, 200))),
                    'bit_score': float(np.random.uniform(300, 600)),
                    'description': vf['description'],
                    'confidence': vf['confidence']
                })
            
            # Additional MEDIUM confidence matches (50% chance per sequence)
            if np.random.random() > 0.5:
                med_genes = ['virulence', 'pathogen', 'toxin', 'secretion']
                selected_gene = np.random.choice(med_genes)
                
                blast_results.append({
                    'query_id': seq_id,
                    'subject_id': f"{selected_gene}_related_gene",
                    'identity': float(np.random.uniform(80, 89)),
                    'alignment_length': int(seq['length'] * np.random.uniform(0.7, 0.85)),
                    'e_value': float(10 ** (-np.random.uniform(50, 99))),
                    'bit_score': float(np.random.uniform(200, 400)),
                    'description': f'{selected_gene.capitalize()} related factor',
                    'confidence': 'Medium'
                })
            
            # LOW confidence hits (reduced to 20% chance)
            if np.random.random() > 0.8:
                blast_results.append({
                    'query_id': seq_id,
                    'subject_id': f"hypothetical_protein_{np.random.randint(1000, 9999)}",
                    'identity': float(np.random.uniform(70, 79)),
                    'alignment_length': int(seq['length'] * np.random.uniform(0.5, 0.7)),
                    'e_value': float(10 ** (-np.random.uniform(20, 49))),
                    'bit_score': float(np.random.uniform(100, 250)),
                    'description': 'Hypothetical protein',
                    'confidence': 'Low'
                })
        
        return sorted(blast_results, key=lambda x: float(x['bit_score']), reverse=True)

    def run_real_blast(self, genome_file, blast_db_path):
        """Run actual BLAST analysis"""
        # Your existing implementation
        with NamedTemporaryFile(delete=False, suffix=".fna") as tmp:
            tmp.write(genome_file.getvalue())
            tmp.flush()
            out_file = tmp.name + "_blast_results.txt"

            cmd = [
                BLASTN_PATH,
                "-query", tmp.name,
                "-db", blast_db_path,
                "-out", out_file,
                "-evalue", "1e-10",
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
            ]

            try:
                subprocess.run(cmd, check=True)
                df = pd.read_csv(out_file, sep="\t", header=None)
                df.columns = [
                    "QueryID", "SubjectID", "Identity(%)", "Align_Length",
                    "Mismatches", "Gap_Openings", "Q_Start", "Q_End",
                    "S_Start", "S_End", "E-value", "BitScore"
                ]
                return df
            except subprocess.CalledProcessError as e:
                st.error(f"❌ BLAST execution failed:\n{e}")
                return None
            finally:
                # Clean up temporary files
                os.remove(tmp.name)
                if os.path.exists(out_file):
                    os.remove(out_file)

def main():
    # Initialize analyzer
    analyzer = StenotrophomonasAnalyzer()
    
    # Header
    st.markdown("""
    <div class="main-header">
        <h1>🧬 Stenotrophomonas maltophilia Virulence Factor DataBase</h1>
        <p>Comprehensive BLAST analysis and virulence gene identification</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Sidebar
    with st.sidebar:
        st.header("📁 File Upload")
        
        # File uploaders
        uploaded_files = st.file_uploader(
            "Upload FASTA files",
            type=['fasta', 'fa', 'fna', 'txt'],
            accept_multiple_files=True,
            help="Upload your genome sequence files in FASTA format"
        )
        
        reference_file = st.file_uploader(
            "Upload Reference Database (Optional)",
            type=['fasta', 'fa', 'fna', 'txt'],
            help="Upload a reference database for comparison"
        )
        
        st.header("⚙️ Analysis Settings")
        
        # Analysis parameters
        min_identity = st.slider("Minimum Identity (%)", 70, 99, 85)
        max_evalue = st.selectbox("Maximum E-value", [1e-10, 1e-20, 1e-50, 1e-100], index=2)
        show_hypothetical = st.checkbox("Show hypothetical proteins", True)
        use_real_blast = st.checkbox("Use Real BLAST (requires configured BLASTN_PATH)", False)
        
        # Sample data option
        if st.button("Load Sample Data"):
            st.session_state['sample_data'] = True

    # Main content
    if uploaded_files or st.session_state.get('sample_data', False):
        # Process uploaded files or sample data
        if st.session_state.get('sample_data', False):
            # Sample data
                        # Sample data for conversion
            sample_sequences = [
                {
                    'header': "NZ_LT906480.1:156325-157356 Stenotrophomonas maltophilia strain NCTC10257 chromosome 1, complete sequence",
                    'sequence': "TCAGTGCTGTTCGCGTCTCTCGTACCAGGCTTGTGTACGGTTCACGATCCGCACCACGGACAGCATGACCGGCACCTCGATCAGCACACCAACCACCGTCGCCAGCGCCGCCCCGAATGCACGCCGAACAGTCCTACTGCCGTGGCAACCGCCAGCTCGAAAAAATTGCTGGCGCCAATCAGCGCTGACCGGCGGGCCACGCAATGCGCCACGCCCAGGCGCCGTTGAGCAGTAGGCCAGGCCTGAAGTGAAGATACCTGGATCAGGATCGGCACTGCAATCAGTGCAATCACCAGCGGTTGCTTCACGATCTGCTGGCCCTGGAACCGAACAGCACCACCAGCGTCAGCAGCAGCGCGGACAGGGACAGCGGCCCCAGCCTTCGAAGCGCACGCTGCAGGCCCGCCTCGCCGCTGCCCGAAAGATCCAGTGCCGCAGCAGGGCGGCCACCAGCACCGGCACCACAATGTAGAGACCCACCGACAACAGCAGCGTGTCCCACGGCACAGTGATCGCTGACAGGCCCAGCAGCAGCGCCACCACCGGTGCGTATGCGACCACCATGATCGCGTCGTTCAACGCCACCTGGCTCAGGGTGAAGTTGGCATCGCCGCGGCACAGGTTGCTCCAGACGAACACCATGGCAGTGCAGGGCGCGGCCGCCAGCAGGATCAAACCGCCGATGTAGCTGTCCATCTGTGCTGCGGCAGCCAGTCCGCGAAAGCATGGCGCAGGAACAGCCAGCCCAGCAGCGCCATCGAGAACGGCTTGACCGCCCAGTTGATGAACAGTGTTACCCCGATGCCCTTCCAGTGTTGCTGCAACTGGCGCATCGCGCGGAAGTCGACCTTCAGCAGCATCGGGATGAATCATCAGCCAGATCAGCCGCGCCAGCGCGCGTTGACCTTGGCGACTTCGGCGGACGCCAGCACCGCGAAGCGCCAGGTGGCCAGCAGCGTGCCATCGCAATGCACAGCACAACCCCACAGGGTCCAGGTGCCGCTCGAACAGGCTCAT",
                    'length': 1032,
                    'gc_content': 65.12,
                    'virulence_factors': []
                },
                {
                    'header': "NZ_LT906480.1:525894-526481 Stenotrophomonas maltophilia strain NCTC10257 chromosome 1, complete sequence",
                    'sequence': "TCACCCGGCCGACTTCAAGCCACGCAACAGCATTTCCAGCGCCGTCGTCGCCTGCCGCAGGCGCGCGTCGCCATCCTCGCCCTCTGCAATCCAGAACGCGGCCTCGGCCAGGCTGCCGTAGATCAGCGATGCCAGCGCGTGCGCATCCACCGGCATCGCTACGCCCTGCGCGATCAACTGTTCGATCAGGCGCTGCATCGAGTGCACGCAGTGGCGCTGCGATTCCGGTGAGGCGCCACCCAGCACCGCGCGCGCGTCACGCAGCACGATGCGCTGGATCTCCGGTTCCAGCGCCATCTCCAGATAGGCACGGCAGCGCTGCACGAAACCGTCCCAGGCATCATCCGCGCCATCGCTGATCGCCTGCAGGCGCAGGTCGGTCTCGGCGTCGAGCTGGGCGACCACCGCGGCCAGCAGGCCCTTCTTGTCGCCAAAGTGGTGGTAGAGCGCACCCCGGGTCAGGCCGGCTTGGGCGGTCAGGTCATCCATCGAGGTGGCGGCATAGCCATGCTCGCTGAACACGCGGCGGGCGGTGGCCAGCAGGCTGGCACGGGTTTCTTCCATTTCGGCGCGGGTGCGACGGACCAT",
                    'length': 588,
                    'gc_content': 64.80,
                    'virulence_factors': []
                },
                {
                    'header': "NZ_CP147720.1:3829003-3830394 Stenotrophomonas maltophilia strain CGMCC 1.1788 chromosome, complete genome",
                    'sequence': "ATGCTAGCAGACCCAACCCCCTCCCCGCTGGCCCAGCACTGGCGGCAACGCGACCTGCAGGTGCTGTGGCACCCGTGCACGCAGATGCGTGAACACCCGGACACCCTGCCGCTGGTGCCGATCGCCCGCGGCGAAGGTGCCTGGCTGATCGACCACGATGGCAATCGCTACCTGGATGCGGTCAGCAGCTGGTGGACCAACCTGTTCGGCCATGCCGAACCGCGCATCGGTGGCGCCATCGCCGCCCAGGCCACGCAGCTGGAACAGGTGATGCTGGCCGGCTTCGGCCATGAGCCAGCCATCACCCTGGCCGAACGCCTGCTGGCGCTGGCCCCGCGCCAGCCTGGCCGCGAACCGCTGGCCAAGGTGTTCTACGCCGACAACGGCTCCGCCGGCGTGGAAGTGGCGCTGAAGATGACC TTCCAGTATTTCCAGAACCGTGGGGAACCGCGGCGCACGCGTTTCATCGCACTGGAGAACGGCTACCACGGCGAGACCCTGGGCGCGCTGGCGCTGGGTGACATCCCGCTGTATCGCCGCGTCTATGCGCCGCTGCTGGCCGAAGGGCTGTTCGCGCCGTCGCCTGACGCTTACCTGGCCGAGCCCGGCCAGAGCGCCGCGGATCGCGCGCGCCAGGCCGCCGATGGCCTGGCTACGCTGTTCGACCAGCATCCGGGCGAGATCTGTGCGGTGATCCTGGAGCCACGCCTGCAGTGCGCGGGCGGCATGCGCATGCATGACCCGGTGTACCTGCAGCGGGTGCGCGAACTGTGTGATGCGCACGGTGCATTCATGATCGCCGACGAGATCGCCACCGGTTTCGGCCGCACCGGCACCCTGTTCGCCTGCGAGCAGGCCGGAGTGATGCCAGACCTGATGTGCCTCTCGAAAGGCCTGACCGGCGGCTTCCTGCCGCTGGCCGCCGTGCTGGCGACGCAGGCCCTGTACGACGCGTTCCTCGACGACTCGCGCGAGCGCGCGTTCCAGCACTCGCACAGCTACACCGGTAACCCGCTGGCCTGCGCGGCTGCGCTGGCGACGCTGGACATCTTCCGCGACGACGATGTGATCGCCCGCAACCGTGGCATCGCCTCGGTGATGGGCGCGCTGGCTGCGCCATGAGCGACCATCCGCACGTGGCCGACGTCCGCCAGGCCGGCATGGTGGTGGCCTTCGAGCTGTCGCGGGATGGCAACAAGCGCACGCCGTTCGATCCGGCGCTGCGACTGGGCCTGCACGCCTACAAGGCTGCGCTGAAACGTGGCGTGGTGCTGCGTCCGCTGGGCGACGTGCTGTACTGGATGCCACCCTACTGCGTGGATGACGAACAACTGGAGCTGCTGGCACATACCACGCTGGCGGCGATTGACGAGGCGATTGCATGCGCGTGA",
                    'length': 588,
                    'gc_content': 64.80,
                    'virulence_factors': []
                },
                {
                    'header': "NZ_LT906480.1:1714160-1715341 Stenotrophomonas maltophilia strain NCTC10257 chromosome 1, complete sequence",
                    'sequence': "ATGAGCCAGACCCAAGACACCGCGGCCCCGGCCGCCCCCAACCGCCGCGGCAACCTGCTGCGCGGCCTGTTCGTGATCGTCGTGCTGCTGCTTGCTGCACTGGCGCTGTGGTACTTCATGTTCGGCCGTTGGTTCGAAGAGACCGACGATGCCTACGTGCAGGGCAACCAGGTGCAGATCACCCCGCTGGTGGCCGGTACCGTGGTCGCCATCAACGCCGATGACGGCATGCGCGTGGAGCGCGGCCAGCTGCTGGTGCAGCTGGACCCGTCCGACACCGCGGTGGCACTGCAGCAGGCCGAAGCCAACCTGGCCAAGACGGTACGCCAGACCCGTGGTCTGTACCGCAGCGTGGAAGGCGCCCAGGCGGACTTGAACGCCCGCCAGGTGACCCTGAAGCGCGTGCGCGAAGACTTCGCCCGCCGCAAGGACCTGGCCGCCACCGGCGCCATCTCCAACGAAGAACTGGCCCACGCCCGTGACGAGCTGGCTGCCGCCGAAGCGGCCGTGGCCGGCTCGCGCGAGACCGTCGAGCGCAACCGCGCGCTGGTCGACGACACCGTGATCGCCACCCAGCCGGACGTGCAGGCCGCCGCCGCGCAGCTGCGCCAGGCCTTCCTCAACAACGCCCGCGCCGGCATCGTCGCGCCGGTCACCGGCTACGTCGCCCGTCGTTCGGTGCAGGTCGGCCAGCGCGTGCAGCCGGGCAATGCCCTGATGGCCGTGGTGCCGACCGAGCAGATGTGGGTCGAGGCCAACTTCAAGGAAACCCAGCTGCGCCACATGCGCCTGGGCCAGGAAGTGGAGCTGAAGTCGGACCTGTACGCCGGTGACGTGAAGTACAAGGGCCGCATCCAGAGCCTGGGCCTGGGCACCGGCTCCGCGTTCTCGCTGCTGCCGGCGCAGAACGCCAGCGGCAACTGGATCAAGATCGTGCAGCGCGTGCCGGTGCGCATCGCCATCGATGCCAAGCAGCTGGCGAACACCCGCTGCGCATCGGCCTGTCGATGAAGGCTGAAGTGAGCCTGCGCGACCAGAAGGGCGAAGTGCTGCCGAGCACTCCGGCCAAGGGCACGGTGTTCGACACCGACGTGTATGCCAAGCAGCTGCATGATGCCGACGAGGTGATCCACACGATCATCCAGGGCAACCTGCCGCAGCAGCAGGCGAAGGTGGGCTGA",
                    'length': 588,
                    'gc_content': 64.80,
                    'virulence_factors': []
                },
                {
                    'header': "NZ_LT906480.1:4527359-4528048 Stenotrophomonas maltophilia strain NCTC10257 chromosome 1, complete sequence",
                    'sequence': "ATGCGCAGAGGGAGCGCCCCTATCGTGTCAACCGTGCGCCTAGCTAGCAGTCCACTGGCCTTGGATATCGCCACCATCGATCGTTTCCTGGCGCACAGCCACCGGCGGCGATATCCCACCCGTACCGACGTCTTCCGACCCGGTGACCCGGCGGGAACCCTGTATTACGTTGTCAGCGGCTCGGTTTCAATCATGGCCGAGGAAGATGACGATCGAGAGCTTGTGCTGGGCTACTTCGGTGCCGGCGAGTTCGTTGGCGAAATGGGCCTGTTCGTCGAATCGGACCGCCGCGAGGTGATCCTGCGGACCCGTACCGCCTGCGAACTGGCCGAAATCAGCTACGAGCGCCTGCACCAGCTGTTCCTCGGCCCGCTGTCGGCCGACGCCCCGCGCCTGCTGTATGCGCTGGGCCAGCAGATCTCCAAGCGCCTGCTCGATACCAGCCGCAAGGCCAGTCGCCTGGCATTCCTGGATGTTACCGATCGGATCGTGCGCACACTGCACGATCTGGCGCAGGAACCAGAAGCGATGAGCCATCCGCAGGGCAGCCAGCTGCGCGTGTCACGCCAGGAACTGGCGCGCCTGGTCGGCTGCTCGCGAGAAATGGCCGGCCGCGTGCTGAAGAAGCTGCAGACCGACGGCCTGCTGCATGCCCGCGGCAAGACCGTGGTGCTGTACGGCACCCGTTGA",
                    'length': 588,  # Length of the sequence
                    'gc_content': 64.80,  # GC content percentage
                    'virulence_factors': []  # List of virulence factors (if any)
                },
                {
                    'header': "NZ_CP147720.1:2345766-2347631 Stenotrophomonas maltophilia strain CGMCC 1.1788 chromosome, complete genome",
                    'sequence': "ATGACTGCTACTGCCACTACCGCCCCCTTGCGCGTTGCGCTGGTTGGTAACCCCAACAGCGGCAAGACCGCACTGTTCAACCAGCTGACCGGCAGCCGGCAGAAGGTTGCCAACTACACTGGCGTCACCGTCGAGCGCAAGGAAGGCCGTCTGCGTGCGCCGTCCGGCCGCGAGTTCGCGGTACTCGATCTTCCGGGCGCCTACAGCCTGCATCCGGCCAGCCTGGACGAGGCGATCACCCGTGACCTGTGCCGGGGCTTCTACCCGGGCGAAGCCGCGCCGGATGTGCTGCTGTGCGTGATCGACGCCACCAACCTGCGCCTGCACCTGCGCTTCGCGCTGGAACTGCGCGAGCTGGGCAAGCCGATGGTGGTGGCGCTGAACATGGTTGATGCGGCCCAGCGCCGTGGCATCCAGGTGGACGTGGCCGCGCTGGAGCGCGAGCTGGGCGTGCCGGTAGTGGAGACCGTGGCGGTGCGCAAGCAGGGCGCCAAGGCCCTGGTCGAACGCCTGGATGCGATGGTCCCGCACCTGGATGCACCGGTGCCGGGCCCGGAAGGTGGCATCGATTACCACGCCAAGGTGCGCGAGATCCTGTCGGTGGCCGTGCGCATGCCGGCGCGGACCGCGAAGATCGACGACGCGCTGGACCGCTGGCTGCTGCACCCGGTGTTCGGCCTGATCAGCCTGGCGGTGGTGATGTTCCTGATCTTCCAGGCCGTGTACGCCTGGGCGACGCCGCTGATGGACGCCATCGAGGCCGGTTTCGCCTGGCTCGGCGCCTTCGTCGGCAGCGTGCTGCCGGAAGGCCCGCTGGCCAGCCTGCTGACCGACGGCATCATCGCCGGTGTCGGTGGCGTGGTGGTGTTCCTGCCGCAGATCCTGATCCTGTTCTTCTTCATCCTGGTGCTGGAGGAATCCGGCTACCTGCCGCGTGCGGCGTTCCTGCTCGACCGCATGATGGCCGCTGCCGGCCTGTCGGTCGCTCGTTCATCCCGCTGCTGTCCAGCTTCGCCTGCGCAGTGCCGGGCATCATGTCCACGCGCAGCATCCAGGACCCGCGTGACCGCCTGGCCACGATCCTGGTGGCGCCGCTGATGACCTGTTCCGCACGCCTGCCGGTGTACGCGCTGCTGATCGGTTCGTTCATTCCGCAGAAGACGGTGTGGGGCGTGTTCAACCAGCAGGGCCTGGTGCTGTTCGGCCTGTACGCGGCCGGCATCCTCAGCGCGCTGGCGATGTCGTGGATCATGAAGAAGTGGCGCCGCGACAAGAGCGAGCATCCGCTGATGCTGGAACTGCCGTCGTACCGCCTGCCGCACGTGCGTGACCTGGCGGTCGGCCTGTACGAGCGCGGCATGATCTTCCTCAAGCGCGTTGGCGGCATCATCCTGGCGCTGACCATCCTGCTGTGGGTGCTGCTGTCGTTCCCGGCAGCGCCGGCGGGGGCCACCATGCCCGCGATCGATTACAGCTACGCCGGCCAGATCGGCCATGCGATGGCGGCATTCTTCGCGCCGCTGGGCTTCAACTGGCAGATCTGCATCGCGCTGATCCCGGGCCTGGCCGCGCGTGAAGTGGCGGTGTCTTCACTGGCCACCGTGTATGCTGTCGGCAGCCGATGACGATGCCGCCAGCCAAGCACTGACCCCGCTGATCAGCGATGGCTGGTCGCTGGCCACCGCACTGTCGCTGCTGGTCTGGTACATCTACGCCCCGATGTGCATCTCGACGCTGGCCACCATCAAGCGCGAGACCAACTCGTGGAAGCAGATGGGCTTCGCTGCGTTCTACCTGTTCGCGGCGGCCTACGTGGC",
                    'length': 588,  # Length of the sequence
                    'gc_content': 64.80,  # GC content percentage
                    'virulence_factors': []  # List of virulence factors (if any)
                },
                {
                    'header': "NZ_LT906480.1:1143431-1144525 Stenotrophomonas maltophilia strain NCTC10257 chromosome 1, complete sequence",
                    'sequence': "ATGGCGCACCAGCGCGCCTCGGACCTGTTCATCACTGCGGGCATGCCGCCGGCGATGAAGGTCAACGGCAAGATCTCGCCGATCACGCAGACCCCGCTCACGCCGCAGCAGAGCCGCGACCTGGTCCTCAACGTGATGACCCCGGCGCAGCGCGAGGAATTCGAGAAGACCCACGAGTGCAACTTCGCCATCGGCCTGTCCGGTGTCGGCCGCTTCCGCGTCAGCTGCTTCTACCAGCGCAACCAGGTCGGCATGGTGCTGCGTCGCATCGAGACGCGCATTCGACGGTGGAAGAACTGAGCCTGCCGCCGATCATCAAGACGCTGGCGATGACCAAGCGCGGCATCATCCTGTTCGTCGGTGCCACCGGTACCGGTAAATCGACCTCGCTGGCGGCGATGATCGGTTACCGCAACCAGAACTCGACCGGCCACATCATCACCATCGAAGATCCGATCGAGTTCGTGCACAAGCACGAAGGCTGCATCATCACCCAGCGTGAGGTCGGCATCGATACCGACAGCTGGGAAGCCGCGTTGAAGAACACCCTGCGCCAGGCGCCGGACGTGATCATGATCGGTGAGGTGCGTACCCGCGAAGGCATGGACCACGCCATCGCGTTCGCCGAGACCGGTCACCTGGTGCTGTGCACCCTGCATGCCAACAACGCCAACCAGGCGATGGACCGCATCGTCAACTTCTTCCCGGAAGACCGCCGCAACCAGCTGCTGATGGACCTGTCGCTGAACCTCAAGGGCGTGGTCGCGCAGCAGCTGGTGCCGTCGCCCGATGGCCGCTCGCGAAAGGTCGCCATGGAGATCCTGCTCGGCACGCCGCTGGTTCAGGACTACATCCGCGACGGCGAGATCCACAAGCTGAAGGAAGTGATGAAGGACTCGGTCCAGCTGGGCATGAAGACCTTCGACCAGAGCCTGTTCGAGCTGTACCAGGCCGGCGAGATCAGCTACGAGGACGCGCTGCGTTACGCCGATTCGCAGAACGAAGTGCGCCTGCGCATCAAGCTCAGCCAGGGCGGCGACGCACGCACGCTGTCACAGGGCCTGGATGGCGTGGAGATCTCCGAGATCCGGTGA",
                    'length': 588,  # Length of the sequence
                    'gc_content': 64.80,  # GC content percentage (calculated)
                    'virulence_factors': []  # List of virulence factors (if any)
                },
                {
                    'header': "NZ_LT906480.1:3121727-3122146 Stenotrophomonas maltophilia strain NCTC10257 chromosome 1, complete sequence",
                    'sequence': "TTATGCACTGTAGTACTGGACGGCGTGCTTGATCATTTCAAACGCCATGCCAGCGCCGTCGCCACCACCGTTATGCACTGTAGTACTGGACGGCGTGCTTGATCATTTCAAACGCCATGCCAGCGCCGTCGCCACCACCGACCTTTTCGCAGCTGCCAGCCTCCTTGCCGAGCAGCATGCGCGCGATCTGGCCTGCCGCGCCGATCAGCACACCGCCGATGAGGATCGGGGCCACGTCACCGATGCGCTTGTGGGCAAACGCAATCTGGTAACCGGCGAAGATCACGGCGATCGTGACCACGGCGATCGAGGCCATGTTCAGCAGACCGTTGATGTTGGAGAAGAAGCCGCAGACCTTGCCATCGGTACCACCGAAGTCAGTGCCGGCCAGGGCCTGGGGAGCGAACACAGCACCGGCGAATGCGACGGCCATCAGCATGGTCTTCAGCGTGCGCTGGGCCTGGACGAGGTCGAGATTGAATCGCTTCAT",
                    'length': 420,  # Length of the sequence
                    'gc_content': 61.90,  # GC content percentage (calculated)
                    'virulence_factors': []  # List of virulence factors (if any)
                },
                {
                    'header': "NZ_LT906480.1:1342921-1343889 Stenotrophomonas maltophilia strain NCTC10257 chromosome 1, complete sequence",
                    'sequence': "TCAGTACGTGCGGCGCAGTTCCAGCGCGCCATCCGGGCGCGTCTGGCCAGCAATGGACCAGCCAAGTTCCAGCGCGGACAGCAACTCCTGTGCGTCGCCTGCGCGGAACGTGCCGCTGACTGCCAGATCGGTGATCTCCGGGTCGGCGATCACCAGCGGTGCATGGCCATAGCGGTTCATGCGTTCGACCACGATCGACAGCGGGGTGGCATCGAAGACCAGCTGCCCATGCAGCCAGCCCTCGGCGGCTGCCGTGGAGGACAGCGCTGGACCGGGCTGGATGCGCCCGGAGGGGAGTACCTGCAGTTGCTGGCCAGGGGCCAGGGTGTGCTGGGCGGCGCCGTTGCTGACTTCCACCGCACCCTCCAGCAGCGCGACCTCGACCAGGCCATCGTGCAGGCGCTCGACCTGGAAGGTGGTGCCGATATCGCGGATGGTGCTGGTGCCGGCGCGCAGCTGCAGTGCTTTGTTCGATGGCGCGACCTGCAACTGCAGGCGGCCGCGTTCGAGGTCCAGCTCGCGGTGGCGCCAGCCAAAGCGTGCGTGCAGCGTCGTGCCGGCATCGAGCATGGCTACGCTGCCATCGGGCAGGGTCAGCCGCTGCGCCTGATGCGTGTTATTGGCAAAGCGCTGCGGGGCTGGGTTGCCATCCACCGCTACCATCCAGCCCACGCCGATGGCCAGGCAGATACCAGCGGCCGCGGCCGCAGCAGGCAGCCAGCGGCGTCGGCCAGGCCGTGCCGCACGCGCCGCAGCCGTCCGCAGCCATGGGTCCGCAGCCAGGCCCTGGCCCTGCTGGAAAAGATCTTCGGCCTGCGCCCAGGCTTTGACATGGGCAGGGTCTTCAGCCAGCCAGTCTTCGAACGACTCGCGTTCGGCAGGGGTGCAGTCCGGGGCCTCCAACCGCGCCACCCAGGCGCTGGCGCGTGCGAACAGCGAGGCGTCTTCCGTGTGCCGGTGGCTGCTCATTCTTCCGTGTGCCGGTGGCTGCTCAT",
                    'length': 588,  # Length of the sequence
                    'gc_content': 61.90,  # GC content percentage (calculated)
                    'virulence_factors': []  # List of virulence factors (if any)
                },
                {
                    'header': "NZ_CP147720.1:1086924-1087868 Stenotrophomonas maltophilia strain CGMCC 1.1788 chromosome, complete genome",
                    'sequence': "GTGACCGTATCCACCCTCGCTTTTGTCTTCCCCGGCCAGGGCTCGCAGTCTGTGGGCATGGTGGCCGAGCTGGCCGAACTGCACCCGCAGGTGCGTGAAGCCTTCACCGAAGCATCCGATGGTGCCGGCGTCGACCTGTGGGCGCTGTCCCAGGGTGGCCCGGAGGAAATGCTCAACCGTACCGAATACACCCAGCCGGCCCTGCTGGCCGCGAGCATCGGCGTGTGGCGCGCCTGGAACGCCGTGGGCGGTCCGCGCCCGTCGGTGCTGGCCGGCCACAGCCTGGGCGAGTACACCGCCCTGGTCGCCGCTGGCGCGCTGAGCCTGCATGACGGTGCGCACCTGGTGCGCCTGCGCGGCCAGCTGATGCAGGAAGCCGCACCGGCCGGTGTCGGCGCCATGGCCGCCGTGCTCGGTGCTGAAGACCAGCTGGTGCTGGACGTCTGCGCCGAAGCCGCGGGCAGCCAGGTGGTGGTGCCGGCCAACTTCAATTCGCCGGGCCAGATCGTGATCGGCGGCGACGCAGATGCGGTTGACCGCGCGCTGGCGCTGCTGGCCGGGAAGGGCGTGCGCAAGGCCGTCAAGCTGGCCGTGAGCGTGCCCTCGCACACCCCGCTGATGCGCGAAGCCGCCAACCGCCTCGCCGAAGTGATGGCCGGCCTGTCCTGGCAGGCGCCGCAGCTGCCGGTGGTGCAGAACGTGGATGCCCAGGTGCATGACGGTATCGACGCGATCCGTACCGCGCTGGTCCAGCAGCTGTACCAGCCGGTGCAGTGGACCGGTTGCGTGCAGGCGCTGGCCGCCCGTGGCATCACCCAGGTTGCCGAGTGCGGCCCGGGC AAGGTGCTGACCGGCCTGGTCAAGCGCATCGACAAGGCCATCGACGGCCGTTCACTGGCCACCCCGGGCGACTTCGAAGCCGCCCGCGAGGCATGGTCGGCCTGA",
                    'length': 588,  # Length of the sequence
                    'gc_content': 61.90,  # GC content percentage (calculated)
                    'virulence_factors': []  # List of virulence factors (if any)
                },
                {
                    'header': "NZ_LT906480.1:223102-223710 Stenotrophomonas maltophilia strain NCTC10257 chromosome 1, sme beta-lactamase",
                    'sequence': "ATGGCCTACAAACGCTCTGCACTGATGGAAGAACGCCTGGCCGGGCCGTGAACGGATCCTGCTGGCCACCCGCGAGCTGGTCGCCACCGGCGGCTGGCGCAACGCCCCGTGACCGCCGTGGCAACCCAGCCGGGGTATCCACCGGCCTGATCTATCGCCATTTCCCGTCCAAGGCCGAGCTGTTCGTGGAGGTGCTCAACGCCGCGGTGGCCCACGAGGTGGCGATCATGGAACGCATCGCCAGCGGTCCAGAGGCGGCCAGCGAGCGGCTGCGGCTGGCGATCACCGCCTTCGTACGCCGCGCCCTGGCCGGGCCGGGGCTGGCCCATGCCTTCATCGTCGAGCCGGTCGATCCGGACGTGGAGCCGAACGCATGCGTGGCCGTCGTGCCTTTGGGGATGTCTTCCTGCGGCTGGTGGAGGAGGGCGTGGCCGCCGGTGAGCTGCCGGCGCAGGACGCGCATGGGCCGCGCCCTGCCTGGTCGGCGCCTTCACCGAGGCGATGGTCGGCCCGACCGCACCCAGCCGCGAAGCGCACCGCGACGAAGACGCGCTGGTGGATGCGATCTGCAGCTTCTGCCTGCGCGCGATCGGCGCCCGGTAA",
                    'length': 609,
                    'gc_content': 68.97,
                    'virulence_factors': []
                },
                {
                    'header': "NZ_CP147720.1:2418793-2419536 Stenotrophomonas maltophilia strain CGMCC 1.1788 chromosome, complete genome",
                    'sequence': "TCATTCTTCCAGGCCGGCGTCCGCCGCCTCGAAAATCTTGAGCCGGCCACGCAGGCGCAGCACCGCCTGGCCGTGGATCTGGCAGACCCGCGACTCGCTCACGCCGAGCACTGCGCCGATCTCCTTCAGGTTCAGTTCCTGCTCGTAGTAGAGCGAGAGCACCAGCTGTTCGCGCTCGGGCAGGTGGCCGATGGCCTTGCCCAGCTCGCGGCCGAACTCACCACGCTCCAGTACCTGCTGCGGGGTCGGGCCGCCCTGGGCGACGGTGTCCAGCTCGCCCTGGTCCTCGATGCGCGACTCCAGGCTCAGCACCTGGCCACGTGCGGCGTCTTCCATCAACCGCAGGTATTCGGGCAGCGGCATCTCCATCGCGGCGGCCACTTCGGTGGCGCTGGCGGCGCGGCCGCTGCTCTGTTCCAGGCGGCGGATGGTGGCGGCGGCATCGCGTGCACGGCGGTGCACGGAGCGCGGCACCCAGTCGCCACGGCGGATCTCGTCGATCATCGAACCCCGGATGCGGATCGAGGCATAGGTCTCGAACGAGGCGCCCTGGTCGGCGTCGTAGCTGCGCGACGCTTCGATCAGGCCCATCATGCCGGCCTGGATCAGGTCGTCCACTTCAACGCTGGCCGGCAACCGCGCGGCCAGGTGGTGGGCGATGCGCCGCACCAGGTCCGAATGCTGGGCGATGACCTCGTTGGCCGCCGAGCGCTGGACTTCCCTGTACTGGGCTGCGCCTTTCAT",
                    'length': 588,  # Length of the sequence
                    'gc_content': 61.90,  # GC content percentage (calculated)
                    'virulence_factors': []  # List of virulence factors (if any)
                },
                {
                    'header': "NZ_CP147720.1:3862943-3863491 Stenotrophomonas maltophilia strain CGMCC 1.1788 chromosome, complete genome",
                    'sequence': "ATGGCGAAGACCAAGAGCACGACCAAGGCCAAGACCGGCAAGCAGAAGCTTGCTGCGGCGGCACCGTCGGCACCGAACATCGATATCGGGATCACCCAGGGCGACCGCAAGAAGATCGCCGACGGTCTGTCGCGCTTCCAGGCCGATGCGTTCACGCTCTACCTGAAGACGCACAACTTCCACTGGAACGTGACCGGGTCGATGTTCAACTCGCTGCACACCATGTTCGAGACGCAGTACACCGAGCAGTGGGCGGCACTGGACGACGTGGCCGAGCGCATCCGCGCGCTGGGCTTCAATGCCCCGGGCTCCTACCGGGAATTCGCTGCGCTGACCTCGATCGCCGAGGAACCGGGCCTGACCGACAGTGCGGACTGGCGTGAAATGGTACGCCAGTTGGTGGTCGCCAACGAAGCGGTCTGCCGTACGGCACGTGAAGTGCTGGAGGTGGCGGCCAAGGGGGACGACGCCCCGACCGAGGATCTGATGACCCAGCGGCTGCAGACACACGAGAAGTACGCCTGGATGCTGCGTTCTCTGCTCCAGTAA",
                    'length': 588,  # Length of the sequence
                    'gc_content': 61.90,  # GC content percentage (calculated)
                    'virulence_factors': []  # List of virulence factors (if any)
                }







            ]
            
            # Process sample sequences
            sequences = []
            for seq in sample_sequences:
                processed_seq = seq.copy()
                processed_seq['gc_content'] = analyzer.calculate_gc_content(seq['sequence'])
                processed_seq['virulence_factors'] = analyzer.identify_virulence_factors(seq['header'], seq['sequence'])
                sequences.append(processed_seq)
            
        else:
            # Process uploaded files
            sequences = []
            for uploaded_file in uploaded_files:
                file_sequences = analyzer.parse_fasta(uploaded_file.read())
                sequences.extend(file_sequences)
        
        if sequences:
            # Analysis tabs
            tab1, tab2, tab3, tab4 = st.tabs(["📊 Overview", "🧬 Sequences", "🎯 BLAST Results", "📤 Export"])
            
            with tab1:
                st.header("Analysis Overview")
                
                # Statistics
                total_sequences = len(sequences)
                total_virulence = sum(len(seq['virulence_factors']) for seq in sequences)
                avg_length = int(np.mean([seq['length'] for seq in sequences]))
                avg_gc = round(np.mean([seq['gc_content'] for seq in sequences]), 2)
                
                # Display metrics
                col1, col2, col3, col4 = st.columns(4)
                
                with col1:
                    st.metric("Total Sequences", total_sequences)
                
                with col2:
                    st.metric("Virulence Genes Found", total_virulence)
                
                with col3:
                    st.metric("Average Length (bp)", f"{avg_length:,}")
                
                with col4:
                    st.metric("Average GC Content", f"{avg_gc}%")
                
                # Visualizations
                col1, col2 = st.columns(2)
                
                with col1:
                    # Length distribution
                    lengths = [seq['length'] for seq in sequences]
                    fig = px.histogram(
                        x=lengths,
                        title="Sequence Length Distribution",
                        labels={'x': 'Length (bp)', 'y': 'Count'},
                        nbins=20
                    )
                    fig.update_layout(showlegend=False)
                    st.plotly_chart(fig, use_container_width=True)
                
                with col2:
                    # GC content distribution
                    gc_contents = [seq['gc_content'] for seq in sequences]
                    fig = px.histogram(
                        x=gc_contents,
                        title="GC Content Distribution",
                        labels={'x': 'GC Content (%)', 'y': 'Count'},
                        nbins=20
                    )
                    fig.update_layout(showlegend=False)
                    st.plotly_chart(fig, use_container_width=True)
                
                # Virulence factors summary
                st.subheader("Virulence Factors Summary")
                
                virulence_counts = {}
                for seq in sequences:
                    for vf in seq['virulence_factors']:
                        gene = vf['gene']
                        if gene in virulence_counts:
                            virulence_counts[gene] += 1
                        else:
                            virulence_counts[gene] = 1
                
                if virulence_counts:
                    vf_df = pd.DataFrame(
                        list(virulence_counts.items()),
                        columns=['Gene', 'Count']
                    ).sort_values('Count', ascending=False)
                    
                    fig = px.bar(
                        vf_df,
                        x='Gene',
                        y='Count',
                        title="Virulence Factors Distribution"
                    )
                    st.plotly_chart(fig, use_container_width=True)
                else:
                    st.info("No virulence factors identified in the sequences.")
            
            with tab2:
                st.header("Sequence Analysis")
                
                # Create a DataFrame for sequences
                seq_data = []
                for i, seq in enumerate(sequences):
                    seq_data.append({
                        'Sequence ID': f'Sequence_{i+1}',
                        'Header': seq['header'],
                        'Length (bp)': seq['length'],
                        'GC Content (%)': seq['gc_content'],
                        'Virulence Factors': "; ".join([vf['gene'] for vf in seq['virulence_factors']])
                    })
                
                seq_df = pd.DataFrame(seq_data)
                st.dataframe(seq_df, use_container_width=True)  # Display sequences in a table format

                for i, seq in enumerate(sequences[:5]):  # Show first 5 with expanders
                    with st.expander(f"Sequence Detail: {seq['header'][:100]}...", expanded=(i<2)):
                        st.write(f"**Full Header:** {seq['header']}")
                        st.write(f"**Length:** {seq['length']} bp")
                        st.write(f"**GC Content:** {seq['gc_content']}%")
                        
                        if seq['virulence_factors']:
                            st.write("**Virulence Factors:**")
                            for vf in seq['virulence_factors']:
                                st.markdown(f"""
                                <div class="virulence-factor">
                                    {vf['gene']}: {vf['description']} ({vf['confidence']} confidence)
                                </div>
                                """, unsafe_allow_html=True)
                        
                        st.write("**Sequence Preview (first 200bp):**")
                        st.code(seq['sequence'][:200] + ("..." if len(seq['sequence']) > 200 else ""))
            
            with tab3:
                st.header("BLAST Analysis Results")
                
                if st.button("Run BLAST Analysis", type="primary"):
                    with st.spinner("Running BLAST analysis..."):
                        if use_real_blast:
                            if reference_file is not None:
                                result_df = analyzer.run_real_blast(reference_file, BLAST_DB)
                                if result_df is not None:
                                    st.session_state['blast_results'] = result_df.to_dict('records')
                            else:
                                st.error("No reference file available for real BLAST analysis")
                        else:
                            blast_results = analyzer.simulate_blast_analysis(sequences)
                            st.session_state['blast_results'] = blast_results
                
                if 'blast_results' in st.session_state:
                    blast_results = st.session_state['blast_results']
                    
                    # Create a DataFrame for BLAST results
                    blast_data = []
                    for result in blast_results:
                        row = {
                            'Query ID': result['query_id'],
                            'Subject ID': result['subject_id'],
                            'Identity (%)': f"{result['identity']:.1f}%",
                            'Alignment Length': result.get('alignment_length', 'N/A'),
                            'E-value': f"{result['e_value']:.2e}",
                            'Bit Score': f"{result['bit_score']:.1f}",
                            'Description': result['description'],
                            'Confidence': result.get('confidence', 'N/A')
                        }
                        blast_data.append(row)
                    
                    blast_df = pd.DataFrame(blast_data)
                    st.subheader(f"Found {len(blast_results)} matches")
                    st.dataframe(blast_df, use_container_width=True)  # Display BLAST results in a table format
                    
                    # Visualization for top hits
                    if len(blast_results) > 0:
                        st.subheader("Top Matches Visualization")
                        
                        # Show top 15 hits for better visualization
                        viz_df = blast_df.head(15)
                        fig = px.bar(
                            viz_df,
                            x='Query ID',  # Ensure this matches the DataFrame column name
                            y='Identity (%)',  # Ensure this matches the DataFrame column name
                            color='Description',  # Ensure this matches the DataFrame column name
                            title="Top BLAST Matches by Identity %",
                            labels={'Identity (%)': 'Identity (%)'}
                        )
                        st.plotly_chart(fig, use_container_width=True)
                    
                    if len(blast_results) > 20:
                        st.info(f"Showing {len(blast_results)} matches. Use the table filters to explore.")
                else:
                    st.info("Click 'Run BLAST Analysis' to see results.")

            with tab4:
                st.header("Export Results")
                
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    # CSV Export
                    if st.button("📊 Generate CSV Report"):
                        csv_data = []
                        for i, seq in enumerate(sequences):
                            virulence_list = "; ".join([f"{vf['gene']}:{vf['description']}" for vf in seq['virulence_factors']])
                            csv_data.append({
                                'Sequence_ID': f'Sequence_{i+1}',
                                'Header': seq['header'],
                                'Length_bp': seq['length'],
                                'GC_Content_%': seq['gc_content'],
                                'Virulence_Factors': virulence_list,
                                'Virulence_Count': len(seq['virulence_factors'])
                            })
                        
                        if 'blast_results' in st.session_state:
                            for result in st.session_state['blast_results']:
                                csv_data.append({
                                    'Sequence_ID': result['query_id'],
                                    'Subject_ID': result['subject_id'],
                                    'Identity_%': result['identity'],
                                    'Alignment_Length': result.get('alignment_length', 'N/A'),
                                    'E-value': result['e_value'],
                                    'Bit_Score': result['bit_score'],
                                    'Description': result['description'],
                                    'Confidence': result.get('confidence', 'N/A')
                                })
                        
                        df = pd.DataFrame(csv_data)
                        csv = df.to_csv(index=False)
                        
                        st.download_button(
                            label="Download CSV",
                            data=csv,
                            file_name=f"stenotrophomonas_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                            mime="text/csv"
                        )
                
                with col2:
                    # FASTA Export
                    if st.button("🧬 Generate FASTA File"):
                        fasta_content = ""
                        for i, seq in enumerate(sequences):
                            fasta_content += f">Sequence_{i+1}|{seq['header']}\n"
                            # Format sequence in 80-character lines
                            sequence = seq['sequence']
                            for j in range(0, len(sequence), 80):
                                fasta_content += sequence[j:j+80] + "\n"
                        
                        st.download_button(
                            label="Download FASTA",
                            data=fasta_content,
                            file_name=f"stenotrophomonas_processed_{datetime.now().strftime('%Y%m%d_%H%M%S')}.fasta",
                            mime="text/plain"
                        )
                
                with col3:
                    # JSON Export
                    if st.button("📋 Generate JSON Report"):
                        json_data = {
                            'analysis_info': {
                                'timestamp': datetime.now().isoformat(),
                                'total_sequences': len(sequences),
                                'total_virulence_genes': sum(len(seq['virulence_factors']) for seq in sequences),
                                'average_length': int(np.mean([seq['length'] for seq in sequences])),
                                'average_gc_content': round(np.mean([seq['gc_content'] for seq in sequences]), 2)
                            },
                            'sequences': sequences
                        }
                        
                        if 'blast_results' in st.session_state:
                            json_data['blast_results'] = st.session_state['blast_results']
                        
                        json_str = json.dumps(json_data, indent=2)
                        
                        st.download_button(
                            label="Download JSON",
                            data=json_str,
                            file_name=f"stenotrophomonas_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                            mime="application/json"
                        )
                
                # Summary report
                st.subheader("Analysis Summary")
                
                summary_report = f"""
## Stenotrophomonas maltophilia Analysis Report

**Analysis Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

### Summary Statistics
- **Total Sequences Analyzed:** {len(sequences)}
- **Total Virulence Genes Found:** {sum(len(seq['virulence_factors']) for seq in sequences)}
- **Average Sequence Length:** {int(np.mean([seq['length'] for seq in sequences])):,} bp
- **Average GC Content:** {round(np.mean([seq['gc_content'] for seq in sequences]), 2)}%

### Virulence Factors Detected
"""
                
                virulence_summary = {}
                for seq in sequences:
                    for vf in seq['virulence_factors']:
                        gene = vf['gene']
                        desc = vf['description']
                        if gene not in virulence_summary:
                            virulence_summary[gene] = {'description': desc, 'count': 0}
                        virulence_summary[gene]['count'] += 1
                
                for gene, info in sorted(virulence_summary.items(), key=lambda x: x[1]['count'], reverse=True):
                    summary_report += f"- **{gene}**: {info['description']} (found in {info['count']} sequences)\n"
                
                if 'blast_results' in st.session_state:
                    summary_report += "\n### BLAST Analysis Summary\n"
                    blast_results = st.session_state['blast_results']
                    summary_report += f"- **Total Matches:** {len(blast_results)}"
                    top_matches = sorted(blast_results, key=lambda x: x.get('bit_score', 0), reverse=True)[:5]
                    if top_matches:
                        summary_report += "\n- **Top Matches:**"
                        for match in top_matches:
                            summary_report += (f"\n  - {match['query_id']} → {match['subject_id']} "
                                             f"(Identity: {match['identity']:.1f}%, "
                                             f"Bit Score: {match['bit_score']:.1f})")
                
                if not virulence_summary:
                    summary_report += "- No virulence factors identified\n"
                
                st.markdown(summary_report)
                
                # Download summary report
                st.download_button(
                    label="📄 Download Summary Report",
                    data=summary_report,
                    file_name=f"stenotrophomonas_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md",
                    mime="text/markdown"
                )
        
        else:
            st.error("No sequences found in uploaded files. Please check your FASTA format.")
    
    else:
        # Welcome screen
        st.markdown("""
        ## Welcome to the S. maltophilia Virulence Gene Analyzer
        
        This tool helps you analyze Stenotrophomonas maltophilia genome sequences to identify virulence factors and perform BLAST analysis.
        
        ### Features:
        - 📁 **Multi-file FASTA upload**
        - 🧬 **Sequence analysis** (length, GC content)
        - 🎯 **Virulence gene identification**
        - 🔍 **Simulated and real BLAST analysis**
        - 📊 **Interactive visualizations**
        - 📤 **Multiple export formats**
        
        ### Getting Started:
        1. Upload your FASTA files using the sidebar
        2. Or click "Load Sample Data" to try example sequences
        3. Explore analysis results in the tabs
        4. Export your results
        """)

if __name__ == "__main__":
    main()
