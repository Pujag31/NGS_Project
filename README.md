# VCF GENERATION
import streamlit as st
import subprocess
import os
import zipfile
import pandas as pd

from reportlab.platypus import (
    SimpleDocTemplate,
    Paragraph,
    Spacer
)

from reportlab.lib.styles import (
    getSampleStyleSheet
)

# ---------------------------------

st.set_page_config(
    page_title="NGS QC + Variant Dashboard",
    layout="wide"
)
st.markdown(
"""
<style>

.main{

background:
linear-gradient(
135deg,
#f8fbff,
#edf6ff
);

}

h1{

text-align:center;

color:#0047AB;

font-size:44px;

}

.stMetric{

background:white;

padding:20px;

border-radius:18px;

box-shadow:
0 4px 18px
rgba(0,0,0,0.08);

}

.stButton>button{

width:100%;

border-radius:15px;

height:52px;

font-size:18px;

}

.stDownloadButton>button{

background:
linear-gradient(
90deg,
#1f77ff,
#5f8dff
);

color:white;

border-radius:15px;

height:50px;

}

</style>
""",
unsafe_allow_html=True
)

# ---------------------------------
# FASTQC PARSER
# ---------------------------------

def parse_fastqc(path):

    metrics = {
        "Reads": "-",
        "Length": "-",
        "GC": "-"
    }

    with open(path) as f:

        for line in f:

            line = line.strip()

            if line.startswith(
                "Total Sequences"
            ):

                metrics["Reads"] = (
                    line.split("\t")[1]
                )

            elif line.startswith(
                "Sequence length"
            ):

                metrics["Length"] = (
                    line.split("\t")[1]
                )

            elif line.startswith(
                "%GC"
            ):

                metrics["GC"] = (
                    line.split("\t")[1]
                )

    return metrics


# ---------------------------------
# PDF
# ---------------------------------

def generate_pdf(df):

    pdf = (
        "output/QC_Report.pdf"
    )

    doc = (
        SimpleDocTemplate(
            pdf
        )
    )

    styles = (
        getSampleStyleSheet()
    )

    report = []

    report.append(

        Paragraph(

            "NGS QC REPORT",

            styles[
                "Title"
            ]

        )

    )

    report.append(
        Spacer(
            1,
            20
        )
    )

    for _, row in df.iterrows():

        gc = str(
            row["GC"]
        )

        status = (
            "PASS"
        )

        try:

            gc = int(gc)

            if gc < 40:

                status = (
                    "WARN"
                )

            if gc > 60:

                status = (
                    "WARN"
                )

        except:

            status = (
                "UNKNOWN"
            )

        txt = f"""

        <b>Sample</b><br/>
        {row["Sample"]}

        <br/><br/>

        <b>Total Reads</b><br/>
        {row["Reads"]}

        <br/><br/>

        <b>Sequence Length</b><br/>
        {row["Length"]}

        <br/><br/>

        <b>GC Content</b><br/>
        {row["GC"]} %

        <br/><br/>

        <b>QC Status</b><br/>
        {status}

        """

        report.append(
            Paragraph(
                txt,
                styles[
                    "BodyText"
                ]
            )
        )

        report.append(
            Spacer(
                1,
                15
            )
        )

    doc.build(
        report
    )

    return pdf


# ---------------------------------
st.markdown(
"""
<h1>
🧬 GeneScope
</h1>

<h4 style='text-align:center;'>

Automated NGS QC +
Variant Analysis Dashboard

</h4>

""",

unsafe_allow_html=True
)

# ---------------------------------

if files:

    os.makedirs(
        "uploads",
        exist_ok=True
    )

    os.makedirs(
        "output",
        exist_ok=True
    )

    rows = []

    for file in files:

        save = (
            "uploads/"
            + file.name
        )

        with open(
            save,
            "wb"
        ) as f:

            f.write(
                file.getbuffer()
            )

        subprocess.run(
            [
                "fastqc",
                save,
                "-o",
                "output"
            ]
        )

        zip_file = (
            file.name
            .split(".")[0]

            + "_fastqc.zip"
        )

        zip_path = (
            "output/"
            + zip_file
        )

        if os.path.exists(
            zip_path
        ):

            with zipfile.ZipFile(
                zip_path
            ) as z:

                target = [

                    x

                    for x

                    in z.namelist()

                    if (
                        "fastqc_data.txt"
                        in x
                    )

                ][0]

                z.extract(
                    target,
                    "output"
                )

                data = (
                    parse_fastqc(
                        "output/"
                        + target
                    )
                )

                data[
                    "Sample"
                ] = (
                    file.name
                )

                rows.append(
                    data
                )

    # ----------------
    # DASHBOARD
    # ----------------

    st.header(
        "📊 QC Dashboard"
    )

    df = (
        pd.DataFrame(
            rows
        )
    )

    st.dataframe(
        df,
        use_container_width=True
    )

    c1, c2, c3 = (
        st.columns(3)
    )

    c1.metric(
        "Samples",
        len(df)
    )

    c2.metric(
        "Average Reads",
        df[
            "Reads"
        ].iloc[0]
    )

    c3.metric(
        "Average GC",
        df[
            "GC"
        ].iloc[0]
    )

    # ----------------
    # PDF
    # ----------------

    pdf = (
        generate_pdf(
            df
        )
    )

    st.subheader(
        "⬇ Downloads"
    )

    with open(
        pdf,
        "rb"
    ) as f:

        st.download_button(
            "Download QC PDF",
            f,
            "QC_Report.pdf"
        )

    # ----------------
    # PIPELINE
    # ----------------

    first = (
        files[0]
        .name
    )

    fastq = (
        "uploads/"
        + first
    )

    st.write(
        "Generating VCF..."
    )

    subprocess.run(
        f"""
        bwa mem
        reference/ref.fa
        {fastq}
        > output/aligned.sam
        """,
        shell=True
    )

    subprocess.run(
        """
        samtools view
        -Sb
        output/aligned.sam
        > output/aligned.bam
        """,
        shell=True
    )

    subprocess.run(
        [
            "samtools",
            "sort",
            "output/aligned.bam",
            "-o",
            "output/sorted.bam"
        ]
    )

    subprocess.run(
        [
            "samtools",
            "index",
            "output/sorted.bam"
        ]
    )

    subprocess.run(
        """
        bcftools mpileup
        -f reference/ref.fa
        output/sorted.bam |

        bcftools call
        -mv
        -Ov
        -o output/variants.vcf
        """,
        shell=True
    )

    # SNP

    subprocess.run(
        """
        bcftools view
        -v snps
        output/variants.vcf
        -o output/snps.vcf
        """,
        shell=True
    )

    # INDEL

    subprocess.run(
        """
        bcftools view
        -v indels
        output/variants.vcf
        -o output/indels.vcf
        """,
        shell=True
    )

    # ----------------
    # DOWNLOADS
    # ----------------

    if os.path.exists(
        "output/variants.vcf"
    ):

        with open(
            "output/variants.vcf",
            "rb"
        ) as f:

            st.download_button(
                "Download Full VCF",
                f,
                "variants.vcf"
            )

    if os.path.exists(
        "output/snps.vcf"
    ):

        with open(
            "output/snps.vcf",
            "rb"
        ) as f:

            st.download_button(
                "Download SNP VCF",
                f,
                "snps.vcf"
            )

    if os.path.exists(
        "output/indels.vcf"
    ):

        with open(
            "output/indels.vcf",
            "rb"
        ) as f:

            st.download_button(
                "Download INDEL VCF",
                f,
                "indels.vcf"
            )

    st.success(
        "Pipeline Completed!"
    )
c1,c2,c3,c4=st.columns(4)

c1.metric(
"🧪 Samples",
len(df)
)

c2.metric(
"🧬 SNPs",
"-"
)

c3.metric(
"🧩 INDELs",
"-"
)

c4.metric(
"📄 Reports",
"Ready"
)
