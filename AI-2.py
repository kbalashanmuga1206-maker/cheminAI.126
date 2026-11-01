import streamlit as st
import io
import json

from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, AllChem
from rdkit.Chem import Crippen, Lipinski, rdMolDescriptors

import streamlit.components.v1 as components
from reportlab.pdfgen import canvas


# PAGE SETTINGS
st.set_page_config(
    page_title="Molecule Studio",
    page_icon="🧪",
    layout="wide"
)

st.title("🧪 Molecule Studio")
st.caption("Simple Chemistry Laboratory")

# NAVIGATION
menu = st.sidebar.selectbox(
    "Choose Tool",
    [
        "Molecule Explorer",
        "Compare Molecules",
        "Chemistry Quiz",
        "About"
    ]
)


# MOLECULE ANALYSIS
def analyze(smiles):

    if not smiles or not smiles.strip():
        return None

    try:
        mol = Chem.MolFromSmiles(smiles.strip())

        if mol is None:
            return None

        properties = {
            "Formula": rdMolDescriptors.CalcMolFormula(mol),
            "Molecular Weight": round(Descriptors.MolWt(mol), 2),
            "LogP": round(Crippen.MolLogP(mol), 2),
            "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 2),
            "H-Bond Donors": Lipinski.NumHDonors(mol),
            "H-Bond Acceptors": Lipinski.NumHAcceptors(mol),
            "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
            "Heavy Atoms": mol.GetNumHeavyAtoms()
        }

        return mol, properties

    except Exception:
        return None


# PDF REPORT
def create_pdf(name, smiles, properties):

    data = io.BytesIO()

    pdf = canvas.Canvas(data)

    pdf.setTitle("Molecule Studio Report")

    pdf.setFont("Helvetica-Bold", 20)
    pdf.drawString(50, 790, "Molecule Studio Report")

    pdf.setFont("Helvetica", 12)

    pdf.drawString(50, 750, f"Molecule: {name}")
    pdf.drawString(50, 725, f"SMILES: {smiles}")

    pdf.setFont("Helvetica-Bold", 14)
    pdf.drawString(50, 685, "Molecular Properties")

    pdf.setFont("Helvetica", 11)

    y = 655

    for key, value in properties.items():

        pdf.drawString(
            60,
            y,
            f"{key}: {value}"
        )

        y -= 25

    pdf.save()

    data.seek(0)

    return data.getvalue()


# 3D VISUALIZATION
def show_3d(mol):

    try:

        mol3d = Chem.AddHs(Chem.Mol(mol))

        status = AllChem.EmbedMolecule(
            mol3d,
            randomSeed=7
        )

        if status != 0:
            st.warning("3D structure generation failed.")
            return

        AllChem.MMFFOptimizeMolecule(mol3d)

        molblock = Chem.MolToMolBlock(mol3d)

        block = json.dumps(molblock)

        html = f"""
        <div id="viewer" style="width:100%;height:400px;"></div>

        <script src="https://3Dmol.org/build/3Dmol-min.js"></script>

        <script>

        const element = document.getElementById("viewer");

        if (window.$3Dmol) {{

            const viewer = $3Dmol.createViewer(
                element,
                {{backgroundColor:"white"}}
            );

            viewer.addModel({block}, "mol");

            viewer.setStyle(
                {{}},
                {{stick:{{}}, sphere:{{scale:0.25}}}}
            );

            viewer.zoomTo();

            viewer.render();

        }} else {{

            element.innerHTML =
                "3D viewer could not load. Check internet connection.";

        }}

        </script>
        """

        components.html(
            html,
            height=420,
            scrolling=False
        )

    except Exception as error:

        st.warning(f"3D visualization error: {error}")


# MOLECULE EXPLORER
if menu == "Molecule Explorer":

    st.header("Molecule Explorer")

    name = st.text_input(
        "Molecule Name",
        "Ethanol"
    )

    smiles = st.text_input(
        "Enter SMILES",
        "CCO"
    )

    if st.button("Analyze Molecule", type="primary"):

        result = analyze(smiles)

        if result is None:

            st.error("Invalid SMILES notation. Try again.")

        else:

            mol, properties = result

            st.success("Molecule successfully analyzed!")

            col1, col2 = st.columns(2)

            with col1:

                st.subheader("2D Structure")

                image = Draw.MolToImage(
                    mol,
                    size=(500, 350)
                )

                st.image(
                    image,
                    use_container_width=True
                )

            with col2:

                st.subheader("Molecular Properties")

                for key, value in properties.items():

                    st.write(f"**{key}:** {value}")

            st.divider()

            st.subheader("3D Molecular Structure")

            show_3d(mol)

            st.divider()

            pdf = create_pdf(
                name,
                smiles,
                properties
            )

            st.download_button(
                "📄 Download PDF Report",
                data=pdf,
                file_name="molecule_report.pdf",
                mime="application/pdf"
            )


# MOLECULE COMPARISON
elif menu == "Compare Molecules":

    st.header("Compare Two Molecules")

    col1, col2 = st.columns(2)

    with col1:

        name1 = st.text_input(
            "First Molecule",
            "Ethanol"
        )

        s1 = st.text_input(
            "First SMILES",
            "CCO"
        )

    with col2:

        name2 = st.text_input(
            "Second Molecule",
            "Acetic Acid"
        )

        s2 = st.text_input(
            "Second SMILES",
            "CC(=O)O"
        )

    if st.button("Compare Molecules", type="primary"):

        r1 = analyze(s1)
        r2 = analyze(s2)

        if r1 is None or r2 is None:

            st.error("Please enter two valid SMILES strings.")

        else:

            m1, p1 = r1
            m2, p2 = r2

            col1, col2 = st.columns(2)

            with col1:

                st.subheader(name1)

                st.image(
                    Draw.MolToImage(m1),
                    use_container_width=True
                )

            with col2:

                st.subheader(name2)

                st.image(
                    Draw.MolToImage(m2),
                    use_container_width=True
                )

            st.subheader("Property Comparison")

            comparison = []

            for key in p1:

                comparison.append({
                    "Property": key,
                    name1: p1[key],
                    name2: p2[key]
                })

            st.dataframe(
                comparison,
                use_container_width=True,
                hide_index=True
            )


# CHEMISTRY QUIZ
elif menu == "Chemistry Quiz":

    st.header("Chemistry Knowledge Test")

    q1 = st.radio(
        "1. What is the chemical formula of water?",
        ["H2O", "CO2", "O2"]
    )

    q2 = st.radio(
        "2. What does SMILES represent?",
        [
            "Molecular structure notation",
            "Temperature",
            "Pressure"
        ]
    )

    q3 = st.radio(
        "3. Which element has atomic number 6?",
        ["Oxygen", "Carbon", "Hydrogen"]
    )

    if st.button("Submit Quiz"):

        score = 0

        if q1 == "H2O":
            score += 1

        if q2 == "Molecular structure notation":
            score += 1

        if q3 == "Carbon":
            score += 1

        st.success(
            f"Your Score: {score}/3"
        )


# ABOUT
elif menu == "About":

    st.header("About Molecule Studio")

    st.write("""
    Molecule Studio is a simple chemistry learning application.

    Features:

    - SMILES validation
    - 2D molecular visualization
    - 3D molecular visualization
    - Molecular property calculation
    - Molecule comparison
    - PDF report generation
    - Chemistry quiz

    Built using Python, Streamlit and RDKit.
    """)

    st.info(
        "Educational chemistry application."
    )


# FOOTER
st.sidebar.markdown("---")

st.sidebar.caption(
    "Molecule Studio | Python + RDKit"
)
