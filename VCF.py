import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import io
import requests
import json
import time
from collections import defaultdict, Counter
import re
from typing import Dict, List, Tuple, Optional, Union
import concurrent.futures
from functools import lru_cache
import gzip
import asyncio
import aiohttp
from datetime import datetime
import logging
from dataclasses import dataclass, field
import tempfile
import os

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Configure Streamlit page
st.set_page_config(
    page_title="Advanced VCF Analysis Suite v3.0 - Real Annotation Edition",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

@dataclass
class VCFMetrics:
    """Data class to store VCF metrics for better type safety"""
    total_variants: int = 0
    chromosomes: int = 0
    chromosome_list: List[str] = field(default_factory=list)
    snps: int = 0
    indels: int = 0
    insertions: int = 0
    deletions: int = 0
    transitions: int = 0
    transversions: int = 0
    ts_tv_ratio: float = 0.0
    qual_mean: float = 0.0
    qual_median: float = 0.0
    qual_std: float = 0.0

@dataclass
class AnnotationResult:
    """Data class for annotation results"""
    variant_id: str
    gene_symbol: str = "Unknown"
    gene_id: str = "Unknown"
    consequence: str = "Unknown"
    impact: str = "Unknown"
    clinical_significance: str = "Unknown"
    allele_frequency: float = 0.0
    allele_frequency_source: str = "Unknown"
    protein_change: str = "Unknown"
    transcript_id: str = "Unknown"
    exon: str = "Unknown"
    intron: str = "Unknown"
    cadd_score: float = 0.0
    sift_score: float = 0.0
    polyphen_score: float = 0.0
    annotation_source: str = "Unknown"
    raw_response: dict = field(default_factory=dict)

class SampleDataGenerator:
    """Generate sample VCF data for testing"""
    
    @staticmethod
    def generate_sample_vcf() -> str:
        """Generate a realistic sample VCF file"""
        header = """##fileformat=VCFv4.2
##fileDate=20231201
##source=SampleGenerator_v1.0
##reference=hg38
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr7,length=159138663>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chrX,length=156040895>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2
"""
        
        # Realistic variant data with known genes
        variants = [
            # BRCA1 variants (chr17)
            "chr17	43094464	rs80357382	C	T	45.2	PASS	DP=87;AF=0.5;AC=1;AN=2	GT:DP:GQ	0/1:87:45	0/0:78:40",
            "chr17	43091434	rs80357906	A	G	52.1	PASS	DP=92;AF=0.5;AC=1;AN=2	GT:DP:GQ	0/1:92:52	0/0:85:48",
            
            # TP53 variants (chr17)
            "chr17	7674220	rs121912651	C	T	48.7	PASS	DP=95;AF=0.5;AC=1;AN=2	GT:DP:GQ	0/1:95:48	0/0:88:42",
            "chr17	7675088	rs121912656	G	A	41.3	PASS	DP=83;AF=0.5;AC=1;AN=2	GT:DP:GQ	0/1:83:41	0/0:76:38",
            
            # APOE variants (chr19) - using chr1 for this example
            "chr1	154426264	rs429358	T	C	55.8	PASS	DP=105;AF=0.25;AC=1;AN=4	GT:DP:GQ	0/1:105:55	0/0:98:50",
            "chr1	154424312	rs7412	C	T	47.9	PASS	DP=89;AF=0.25;AC=1;AN=4	GT:DP:GQ	0/1:89:47	0/0:82:44",
            
            # CFTR variants (chr7)
            "chr7	117307112	rs113993960	C	T	39.4	PASS	DP=76;AF=0.5;AC=1;AN=2	GT:DP:GQ	0/1:76:39	0/0:69:35",
            "chr7	117465784	rs121908769	G	A	43.6	PASS	DP=81;AF=0.5;AC=1;AN=2	GT:DP:GQ	0/1:81:43	0/0:74:40",
            
            # Some indels
            "chr1	155235816	.	ATCG	A	35.2	PASS	DP=72;AF=0.5;AC=1;AN=2	GT:DP:GQ	0/1:72:35	0/0:65:32",
            "chr2	47641559	.	G	GTCA	38.8	PASS	DP=68;AF=0.5;AC=1;AN=2	GT:DP:GQ	0/1:68:38	0/0:61:34",
            
            # Common variants
            "chr1	230845794	rs6025	G	A	61.4	PASS	DP=112;AF=0.1;AC=1;AN=10	GT:DP:GQ	0/1:112:61	0/0:105:58",
            "chr3	46373453	rs1801133	C	T	58.2	PASS	DP=108;AF=0.3;AC=3;AN=10	GT:DP:GQ	0/1:108:58	1/1:101:55",
            
            # X chromosome variants
            "chrX	154411394	rs104894708	G	A	44.1	PASS	DP=85;AF=0.5;AC=1;AN=2	GT:DP:GQ	0/1:85:44	0/0:78:41",
            "chrX	32389644	rs137852512	A	T	36.7	PASS	DP=73;AF=0.5;AC=1;AN=2	GT:DP:GQ	0/1:73:36	0/0:66:33",
            
            # Low quality variants
            "chr1	12345678	.	A	C	15.3	q30	DP=25;AF=0.2;AC=1;AN=5	GT:DP:GQ	0/1:25:15	0/0:23:12",
            "chr2	87654321	.	T	G	18.7	q30	DP=28;AF=0.3;AC=1;AN=3	GT:DP:GQ	0/1:28:18	0/0:25:15",
        ]
        
        return header + "\n".join(variants)

class MyVariantAnnotator:
    """MyVariant.info API integration for variant annotation"""
    
    def __init__(self):
        self.base_url = "https://myvariant.info/v1"
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'VCF-Analysis-Suite/3.0'
        })
        
    def format_variant_hgvs(self, chrom: str, pos: int, ref: str, alt: str) -> str:
        """Format variant in HGVS nomenclature for MyVariant.info"""
        # Remove 'chr' prefix if present
        chrom_clean = chrom.replace('chr', '')
        
        # Handle different variant types
        if len(ref) == 1 and len(alt) == 1:
            # SNV
            return f"{chrom_clean}:g.{pos}{ref}>{alt}"
        elif len(ref) > len(alt):
            # Deletion
            if len(alt) == 1:
                # Simple deletion
                start_pos = pos + 1
                end_pos = pos + len(ref) - 1
                if start_pos == end_pos:
                    return f"{chrom_clean}:g.{start_pos}del"
                else:
                    return f"{chrom_clean}:g.{start_pos}_{end_pos}del"
            else:
                # Complex deletion
                return f"{chrom_clean}:g.{pos}_{pos + len(ref) - 1}delins{alt[1:]}"
        elif len(alt) > len(ref):
            # Insertion
            if len(ref) == 1:
                # Simple insertion
                return f"{chrom_clean}:g.{pos}_{pos + 1}ins{alt[1:]}"
            else:
                # Complex insertion
                return f"{chrom_clean}:g.{pos}_{pos + len(ref) - 1}delins{alt}"
        else:
            # Complex variant
            return f"{chrom_clean}:g.{pos}_{pos + len(ref) - 1}delins{alt}"
    
    def annotate_variant(self, chrom: str, pos: int, ref: str, alt: str) -> AnnotationResult:
        """Annotate a single variant using MyVariant.info"""
        try:
            hgvs_id = self.format_variant_hgvs(chrom, pos, ref, alt)
            variant_id = f"{chrom}:{pos}:{ref}>{alt}"
            
            # Query MyVariant.info
            url = f"{self.base_url}/variant/{hgvs_id}"
            params = {
                'fields': 'cadd,clinvar,dbnsfp,gnomad_exome,gnomad_genome,dbsnp',
                'assembly': 'hg38'
            }
            
            response = self.session.get(url, params=params, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                return self._parse_myvariant_response(variant_id, data)
            else:
                logger.warning(f"MyVariant.info API error for {hgvs_id}: {response.status_code}")
                return AnnotationResult(
                    variant_id=variant_id,
                    annotation_source="MyVariant.info_Error"
                )
                
        except Exception as e:
            logger.error(f"MyVariant annotation error for {variant_id}: {e}")
            return AnnotationResult(
                variant_id=variant_id,
                annotation_source="MyVariant.info_Error"
            )
    
    def _parse_myvariant_response(self, variant_id: str, data: dict) -> AnnotationResult:
        """Parse MyVariant.info response into AnnotationResult"""
        result = AnnotationResult(
            variant_id=variant_id,
            annotation_source="MyVariant.info",
            raw_response=data
        )
        
        try:
            # Extract ClinVar information
            if 'clinvar' in data:
                clinvar = data['clinvar']
                if isinstance(clinvar, dict):
                    if 'clinical_significance' in clinvar:
                        result.clinical_significance = clinvar['clinical_significance']
                    if 'gene' in clinvar and 'symbol' in clinvar['gene']:
                        result.gene_symbol = clinvar['gene']['symbol']
                
            # Extract gnomAD frequencies
            gnomad_af = 0.0
            af_source = "Unknown"
            
            if 'gnomad_exome' in data and 'af' in data['gnomad_exome']:
                gnomad_af = data['gnomad_exome']['af']['af']
                af_source = "gnomAD_exome"
            elif 'gnomad_genome' in data and 'af' in data['gnomad_genome']:
                gnomad_af = data['gnomad_genome']['af']['af']
                af_source = "gnomAD_genome"
            
            result.allele_frequency = gnomad_af
            result.allele_frequency_source = af_source
            
            # Extract CADD scores
            if 'cadd' in data:
                cadd = data['cadd']
                if isinstance(cadd, dict) and 'phred' in cadd:
                    result.cadd_score = cadd['phred']
            
            # Extract dbNSFP predictions
            if 'dbnsfp' in data:
                dbnsfp = data['dbnsfp']
                if isinstance(dbnsfp, dict):
                    # SIFT score
                    if 'sift' in dbnsfp and 'score' in dbnsfp['sift']:
                        sift_score = dbnsfp['sift']['score']
                        if isinstance(sift_score, list):
                            result.sift_score = sift_score[0] if sift_score else 0.0
                        else:
                            result.sift_score = sift_score
                    
                    # PolyPhen score
                    if 'polyphen2' in dbnsfp and 'hdiv' in dbnsfp['polyphen2']:
                        polyphen = dbnsfp['polyphen2']['hdiv']
                        if 'score' in polyphen:
                            pp_score = polyphen['score']
                            if isinstance(pp_score, list):
                                result.polyphen_score = pp_score[0] if pp_score else 0.0
                            else:
                                result.polyphen_score = pp_score
            
            # Extract dbSNP information
            if 'dbsnp' in data:
                dbsnp = data['dbsnp']
                if isinstance(dbsnp, dict) and 'gene' in dbsnp:
                    gene_info = dbsnp['gene']
                    if isinstance(gene_info, dict) and 'symbol' in gene_info:
                        result.gene_symbol = gene_info['symbol']
                        
        except Exception as e:
            logger.error(f"Error parsing MyVariant response for {variant_id}: {e}")
        
        return result

class EnsemblAnnotator:
    """Ensembl REST API integration for variant annotation"""
    
    def __init__(self):
        self.base_url = "https://rest.ensembl.org"
        self.session = requests.Session()
        self.session.headers.update({
            'Content-Type': 'application/json',
            'Accept': 'application/json'
        })
    
    def annotate_variant(self, chrom: str, pos: int, ref: str, alt: str) -> AnnotationResult:
        """Annotate variant using Ensembl VEP REST API"""
        try:
            variant_id = f"{chrom}:{pos}:{ref}>{alt}"
            
            # Format for Ensembl VEP
            chrom_clean = chrom.replace('chr', '')
            vep_input = f"{chrom_clean} {pos} . {ref} {alt} . . ."
            
            # Query Ensembl VEP
            url = f"{self.base_url}/vep/human/hgvs"
            
            # Use HGVS format for better results
            if len(ref) == 1 and len(alt) == 1:
                hgvs_notation = f"{chrom_clean}:g.{pos}{ref}>{alt}"
            else:
                hgvs_notation = f"{chrom_clean}:g.{pos}_{pos + len(ref) - 1}delins{alt}"
            
            data = {
                "hgvs_notations": [hgvs_notation]
            }
            
            response = self.session.post(url, json=data, timeout=15)
            
            if response.status_code == 200:
                vep_results = response.json()
                return self._parse_ensembl_response(variant_id, vep_results)
            else:
                logger.warning(f"Ensembl VEP API error for {variant_id}: {response.status_code}")
                return AnnotationResult(
                    variant_id=variant_id,
                    annotation_source="Ensembl_Error"
                )
                
        except Exception as e:
            logger.error(f"Ensembl annotation error for {variant_id}: {e}")
            return AnnotationResult(
                variant_id=variant_id,
                annotation_source="Ensembl_Error"
            )
    
    def _parse_ensembl_response(self, variant_id: str, vep_results: list) -> AnnotationResult:
        """Parse Ensembl VEP response into AnnotationResult"""
        result = AnnotationResult(
            variant_id=variant_id,
            annotation_source="Ensembl_VEP",
            raw_response=vep_results
        )
        
        try:
            if vep_results and len(vep_results) > 0:
                vep_data = vep_results[0]
                
                # Get the most severe consequence
                if 'most_severe_consequence' in vep_data:
                    result.consequence = vep_data['most_severe_consequence']
                
                # Extract transcript consequences
                if 'transcript_consequences' in vep_data:
                    consequences = vep_data['transcript_consequences']
                    if consequences:
                        # Get the first (usually most relevant) consequence
                        first_consequence = consequences[0]
                        
                        # Gene information
                        if 'gene_symbol' in first_consequence:
                            result.gene_symbol = first_consequence['gene_symbol']
                        if 'gene_id' in first_consequence:
                            result.gene_id = first_consequence['gene_id']
                        
                        # Transcript information
                        if 'transcript_id' in first_consequence:
                            result.transcript_id = first_consequence['transcript_id']
                        
                        # Impact
                        if 'impact' in first_consequence:
                            result.impact = first_consequence['impact']
                        
                        # Protein change
                        if 'hgvsp' in first_consequence:
                            result.protein_change = first_consequence['hgvsp']
                        
                        # Exon/Intron information
                        if 'exon' in first_consequence:
                            result.exon = first_consequence['exon']
                        if 'intron' in first_consequence:
                            result.intron = first_consequence['intron']
                
                # Extract colocated variants (for frequency information)
                if 'colocated_variants' in vep_data:
                    for variant in vep_data['colocated_variants']:
                        if 'frequencies' in variant:
                            freqs = variant['frequencies']
                            # Try to get gnomAD frequency
                            if 'gnomad_exomes' in freqs and 'af' in freqs['gnomad_exomes']:
                                result.allele_frequency = freqs['gnomad_exomes']['af']
                                result.allele_frequency_source = "gnomAD_exomes"
                            elif 'gnomad_genomes' in freqs and 'af' in freqs['gnomad_genomes']:
                                result.allele_frequency = freqs['gnomad_genomes']['af']
                                result.allele_frequency_source = "gnomAD_genomes"
                            break

        except Exception as e:
            logger.error(f"Error parsing Ensembl response for {variant_id}: {e}")
        
        return result

class OptimizedVCFParser:
    """Highly optimized VCF parser using vectorized operations"""
    
    def __init__(self):
        self.chunk_size = 10000
        
    @st.cache_data(show_spinner="Parsing VCF file...")
    def parse_vcf_fast(_self, file_content: str, max_variants: int = None) -> Tuple[List[str], List[str], pd.DataFrame]:
        """Optimized VCF parsing with chunking and vectorized operations"""
        lines = file_content.strip().split('\n')
        header_lines = []
        samples = []
        
        variants_data = []
        
        header_end_idx = 0
        for i, line in enumerate(lines):
            if line.startswith('##'):
                header_lines.append(line)
            elif line.startswith('#CHROM'):
                parts = line.split('\t')
                samples = parts[9:] if len(parts) > 9 else []
                header_end_idx = i + 1
                break
        
        variant_lines = lines[header_end_idx:]
        if max_variants:
            variant_lines = variant_lines[:max_variants]
            
        for line in variant_lines:
            if line.strip() and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 8:
                    qual_value = _self._parse_qual(parts[5])
                    
                    variant_data = {
                        'CHROM': parts[0],
                        'POS': int(parts[1]) if parts[1].isdigit() else 0,
                        'ID': parts[2] if parts[2] != '.' else '',
                        'REF': parts[3],
                        'ALT': parts[4],
                        'QUAL': qual_value,
                        'FILTER': parts[6],
                        'INFO': parts[7],
                        'FORMAT': parts[8] if len(parts) > 8 else '',
                        'SAMPLES': parts[9:] if len(parts) > 9 else []
                    }
                    variants_data.append(variant_data)
        
        df = pd.DataFrame(variants_data)
        
        if not df.empty:
            df = _self._optimize_dtypes(df)
            
        return header_lines, samples, df
    
    @staticmethod
    def _parse_qual(qual_str: str) -> float:
        """Fast QUAL parsing with error handling"""
        if qual_str == '.' or not qual_str:
            return 0.0
        try:
            return float(qual_str)
        except (ValueError, TypeError):
            return 0.0
    
    @staticmethod
    def _optimize_dtypes(df: pd.DataFrame) -> pd.DataFrame:
        """Optimize DataFrame data types for memory efficiency"""
        optimizations = {
            'CHROM': 'category',
            'ID': 'string',
            'REF': 'category',
            'ALT': 'category',
            'FILTER': 'category',
            'FORMAT': 'category'
        }
        
        for col, dtype in optimizations.items():
            if col in df.columns:
                try:
                    df[col] = df[col].astype(dtype)
                except Exception:
                    pass
                    
        return df

class FastStatsCalculator:
    """Optimized statistics calculator using vectorized operations"""
    
    @staticmethod
    @st.cache_data
    def calculate_comprehensive_stats(df: pd.DataFrame) -> VCFMetrics:
        """Calculate comprehensive VCF statistics using vectorized operations"""
        if df.empty:
            return VCFMetrics()
        
        ref_lens = df['REF'].str.len()
        alt_lens = df['ALT'].str.len()
        
        snp_mask = (ref_lens == 1) & (alt_lens == 1)
        insertion_mask = ref_lens < alt_lens
        deletion_mask = ref_lens > alt_lens
        
        metrics = VCFMetrics(
            total_variants=len(df),
            chromosomes=df['CHROM'].nunique(),
            chromosome_list=sorted(df['CHROM'].unique()),
            snps=snp_mask.sum(),
            indels=(insertion_mask | deletion_mask).sum(),
            insertions=insertion_mask.sum(),
            deletions=deletion_mask.sum(),
            qual_mean=df['QUAL'].mean() if 'QUAL' in df.columns else 0.0,
            qual_median=df['QUAL'].median() if 'QUAL' in df.columns else 0.0,
            qual_std=df['QUAL'].std() if 'QUAL' in df.columns else 0.0
        )
        
        transitions, transversions = FastStatsCalculator._calculate_ts_tv(df, snp_mask)
        metrics.transitions = transitions
        metrics.transversions = transversions
        metrics.ts_tv_ratio = transitions / transversions if transversions > 0 else 0.0
        
        return metrics
    
    @staticmethod
    def _calculate_ts_tv(df: pd.DataFrame, snp_mask: pd.Series) -> Tuple[int, int]:
        """Calculate transitions and transversions efficiently"""
        if not snp_mask.any():
            return 0, 0
            
        snp_df = df[snp_mask]
        
        ref_alt = snp_df['REF'].str.upper() + snp_df['ALT'].str.upper()
        
        transition_patterns = {'AG', 'GA', 'CT', 'TC'}
        transversion_patterns = {'AC', 'AT', 'CA', 'CG', 'GC', 'GT', 'TA', 'TG'}
        
        transitions = ref_alt.isin(transition_patterns).sum()
        transversions = ref_alt.isin(transversion_patterns).sum()
        
        return transitions, transversions

class AdvancedAnnotationManager:
    """Manages both real and demo annotations"""
    
    def __init__(self):
        self.myvariant = MyVariantAnnotator()
        self.ensembl = EnsemblAnnotator()
        
    def annotate_variants(self, df: pd.DataFrame, annotation_type: str, max_variants: int = 50) -> pd.DataFrame:
        """Annotate variants using specified method"""
        if annotation_type == "demo":
            return self._demo_annotation(df, max_variants)
        elif annotation_type == "myvariant":
            return self._myvariant_annotation(df, max_variants)
        elif annotation_type == "ensembl":
            return self._ensembl_annotation(df, max_variants)
        elif annotation_type == "combined":
            return self._combined_annotation(df, max_variants)
        else:
            return df.copy()
    
    def _demo_annotation(self, df: pd.DataFrame, max_variants: int) -> pd.DataFrame:
        annotated_df = df.head(max_variants).copy()

        # Mock annotations
        consequences = ['missense_variant', 'synonymous_variant', 'stop_gained', 
                        'splice_variant', 'intron_variant', '3_prime_UTR_variant']
        impacts = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
        clinical_sigs = ['Benign', 'Likely benign', 'Uncertain significance', 
                         'Likely pathogenic', 'Pathogenic', 'Not provided']

        np.random.seed(42)

        annotated_df['gene_symbol'] = [f"GENE_{i+1}" for i in range(len(annotated_df))]
        annotated_df['consequence'] = np.random.choice(consequences, len(annotated_df))
        annotated_df['impact'] = np.random.choice(impacts, len(annotated_df))
        annotated_df['clinical_significance'] = np.random.choice(clinical_sigs, len(annotated_df))
        annotated_df['allele_frequency'] = np.random.uniform(0.001, 0.1, len(annotated_df))
        annotated_df['cadd_score'] = np.random.uniform(0.1, 35.0, len(annotated_df))
        annotated_df['annotation_source'] = 'Demo'

        return annotated_df

    def _myvariant_annotation(self, df: pd.DataFrame, max_variants: int) -> pd.DataFrame:
        annotated_df = df.head(max_variants).copy()
        annotations = []

        progress_bar = st.progress(0)
        status_text = st.empty()

        for i, (_, row) in enumerate(annotated_df.iterrows()):
            try:
                status_text.text(f"MyVariant.info annotation: {i+1}/{len(annotated_df)}")
                result = self.myvariant.annotate_variant(row['CHROM'], row['POS'], row['REF'], row['ALT'])
                annotations.append(result)
                time.sleep(0.1)
            except Exception as e:
                logger.error(f"MyVariant annotation error: {e}")
                annotations.append(AnnotationResult(
                    variant_id=f"{row['CHROM']}:{row['POS']}:{row['REF']}>{row['ALT']}",
                    annotation_source="MyVariant_Error"
                ))
            progress_bar.progress((i + 1) / len(annotated_df))

        self._add_annotations_to_df(annotated_df, annotations)
        progress_bar.empty()
        status_text.empty()
        return annotated_df

    def _ensembl_annotation(self, df: pd.DataFrame, max_variants: int) -> pd.DataFrame:
        annotated_df = df.head(max_variants).copy()
        annotations = []

        progress_bar = st.progress(0)
        status_text = st.empty()

        for i, (_, row) in enumerate(annotated_df.iterrows()):
            try:
                status_text.text(f"Ensembl VEP annotation: {i+1}/{len(annotated_df)}")
                result = self.ensembl.annotate_variant(row['CHROM'], row['POS'], row['REF'], row['ALT'])
                annotations.append(result)
                time.sleep(0.2)
            except Exception as e:
                logger.error(f"Ensembl annotation error: {e}")
                annotations.append(AnnotationResult(
                    variant_id=f"{row['CHROM']}:{row['POS']}:{row['REF']}>{row['ALT']}",
                    annotation_source="Ensembl_Error"
                ))
            progress_bar.progress((i + 1) / len(annotated_df))

        self._add_annotations_to_df(annotated_df, annotations)
        progress_bar.empty()
        status_text.empty()
        return annotated_df

    def _combined_annotation(self, df: pd.DataFrame, max_variants: int) -> pd.DataFrame:
        annotated_df = df.head(max_variants).copy()
        combined_annotations = []

        progress_bar = st.progress(0)
        status_text = st.empty()

        for i, (_, row) in enumerate(annotated_df.iterrows()):
            try:
                status_text.text(f"Combined annotation: {i+1}/{len(annotated_df)}")
                myvariant_result = self.myvariant.annotate_variant(row['CHROM'], row['POS'], row['REF'], row['ALT'])
                time.sleep(0.1)
                ensembl_result = self.ensembl.annotate_variant(row['CHROM'], row['POS'], row['REF'], row['ALT'])

                combined_result = AnnotationResult(
                    variant_id=myvariant_result.variant_id,
                    gene_symbol=ensembl_result.gene_symbol if ensembl_result.gene_symbol != "Unknown" else myvariant_result.gene_symbol,
                    gene_id=ensembl_result.gene_id,
                    consequence=ensembl_result.consequence,
                    impact=ensembl_result.impact,
                    clinical_significance=myvariant_result.clinical_significance,
                    allele_frequency=myvariant_result.allele_frequency,
                    allele_frequency_source=myvariant_result.allele_frequency_source,
                    protein_change=ensembl_result.protein_change,
                    transcript_id=ensembl_result.transcript_id,
                    exon=ensembl_result.exon,
                    intron=ensembl_result.intron,
                    cadd_score=myvariant_result.cadd_score,
                    sift_score=myvariant_result.sift_score,
                    polyphen_score=myvariant_result.polyphen_score,
                    annotation_source="Combined_MyVariant_Ensembl"
                )
                combined_annotations.append(combined_result)
                time.sleep(0.2)
            except Exception as e:
                logger.error(f"Combined annotation error: {e}")
                combined_annotations.append(AnnotationResult(
                    variant_id=f"{row['CHROM']}:{row['POS']}:{row['REF']}>{row['ALT']}",
                    annotation_source="Combined_Error"
                ))
            progress_bar.progress((i + 1) / len(annotated_df))

        self._add_annotations_to_df(annotated_df, combined_annotations)
        progress_bar.empty()
        status_text.empty()
        return annotated_df

    def _add_annotations_to_df(self, df: pd.DataFrame, annotations: List[AnnotationResult]):
        """Add annotation results to DataFrame"""
        required_cols = [
            'gene_symbol', 'gene_id', 'consequence', 'impact',
            'clinical_significance', 'allele_frequency', 'allele_frequency_source',
            'protein_change', 'transcript_id', 'exon', 'intron',
            'cadd_score', 'sift_score', 'polyphen_score', 'annotation_source'
        ]
        for col in required_cols:
            if col not in df.columns:
                df[col] = None

        for i, annotation in enumerate(annotations):
            if i < len(df):
                annotation_cols = {
                    'gene_symbol': annotation.gene_symbol,
                    'gene_id': annotation.gene_id,
                    'consequence': annotation.consequence,
                    'impact': annotation.impact,
                    'clinical_significance': annotation.clinical_significance,
                    'allele_frequency': annotation.allele_frequency,
                    'allele_frequency_source': annotation.allele_frequency_source,
                    'protein_change': annotation.protein_change,
                    'transcript_id': annotation.transcript_id,
                    'exon': annotation.exon,
                    'intron': annotation.intron,
                    'cadd_score': annotation.cadd_score,
                    'sift_score': annotation.sift_score,
                    'polyphen_score': annotation.polyphen_score,
                    'annotation_source': annotation.annotation_source
                }
                for col, value in annotation_cols.items():
                    df.at[df.index[i], col] = value



class AdvancedPlotGenerator:
    """Generate advanced interactive plots for VCF analysis"""
    
    @staticmethod
    def create_comprehensive_dashboard(df: pd.DataFrame, metrics: VCFMetrics, 
                                     annotated_df: pd.DataFrame = None) -> dict:
        """Create comprehensive analysis dashboard"""
        plots = {}
        
        # 1. Variant Distribution by Chromosome
        plots['chrom_dist'] = AdvancedPlotGenerator._create_chromosome_distribution(df)
        
        # 2. Quality Score Distribution
        plots['qual_dist'] = AdvancedPlotGenerator._create_quality_distribution(df)
        
        # 3. Variant Type Distribution
        plots['type_dist'] = AdvancedPlotGenerator._create_variant_type_distribution(df)
        
        # 4. Ti/Tv Analysis
        plots['titv'] = AdvancedPlotGenerator._create_titv_analysis(df, metrics)
        
        # 5. Position Density Heatmap
        plots['pos_heatmap'] = AdvancedPlotGenerator._create_position_heatmap(df)
        
        # 6. Advanced Quality Metrics
        plots['qual_metrics'] = AdvancedPlotGenerator._create_quality_metrics_plot(df)
        
        if annotated_df is not None and not annotated_df.empty:
            # 7. Annotation-based plots
            plots['consequence_dist'] = AdvancedPlotGenerator._create_consequence_distribution(annotated_df)
            plots['impact_analysis'] = AdvancedPlotGenerator._create_impact_analysis(annotated_df)
            plots['frequency_analysis'] = AdvancedPlotGenerator._create_frequency_analysis(annotated_df)
            plots['pathogenicity_scores'] = AdvancedPlotGenerator._create_pathogenicity_scores(annotated_df)
        
        return plots
    
    @staticmethod
    def _create_chromosome_distribution(df: pd.DataFrame):
        """Create chromosome distribution plot"""
        chrom_counts = df['CHROM'].value_counts().sort_index()
        
        fig = px.bar(
            x=chrom_counts.index,
            y=chrom_counts.values,
            title="Variant Distribution by Chromosome",
            labels={'x': 'Chromosome', 'y': 'Number of Variants'},
            color=chrom_counts.values,
            color_continuous_scale='viridis'
        )
        
        fig.update_layout(
            xaxis_title="Chromosome",
            yaxis_title="Number of Variants",
            showlegend=False,
            height=500
        )
        
        return fig
    
    @staticmethod
    def _create_quality_distribution(df: pd.DataFrame):
        """Create quality score distribution"""
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Quality Score Distribution', 'Quality Score Box Plot',
                          'Quality vs Position', 'Quality Score Violin Plot'),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"secondary_y": False}, {"secondary_y": False}]]
        )
        
        # Histogram
        fig.add_trace(
            go.Histogram(x=df['QUAL'], nbinsx=50, name='Quality Distribution'),
            row=1, col=1
        )
        
        # Box plot
        fig.add_trace(
            go.Box(y=df['QUAL'], name='Quality Scores'),
            row=1, col=2
        )
        
        # Quality vs Position scatter
        fig.add_trace(
            go.Scatter(x=df['POS'], y=df['QUAL'], mode='markers',
                      name='Quality vs Position', opacity=0.6),
            row=2, col=1
        )
        
        # Violin plot
        fig.add_trace(
            go.Violin(y=df['QUAL'], name='Quality Distribution'),
            row=2, col=2
        )
        
        fig.update_layout(height=800, title_text="Quality Score Analysis")
        return fig
    
    @staticmethod
    def _create_variant_type_distribution(df: pd.DataFrame):
        """Create variant type distribution pie chart"""
        ref_lens = df['REF'].str.len()
        alt_lens = df['ALT'].str.len()
        
        snps = ((ref_lens == 1) & (alt_lens == 1)).sum()
        insertions = (ref_lens < alt_lens).sum()
        deletions = (ref_lens > alt_lens).sum()
        complex_vars = len(df) - snps - insertions - deletions
        
        labels = ['SNPs', 'Insertions', 'Deletions', 'Complex']
        values = [snps, insertions, deletions, complex_vars]
        
        fig = px.pie(
            values=values,
            names=labels,
            title="Variant Type Distribution",
            color_discrete_sequence=px.colors.qualitative.Set3
        )
        
        fig.update_traces(textposition='inside', textinfo='percent+label')
        return fig
    
    @staticmethod
    def _create_titv_analysis(df: pd.DataFrame, metrics: VCFMetrics):
        """Create Ti/Tv ratio analysis"""
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=(f'Ti/Tv Ratio: {metrics.ts_tv_ratio:.2f}', 'Transition vs Transversion'),
            specs=[[{"type": "indicator"}, {"type": "bar"}]]
        )
        
        # Ti/Tv ratio gauge
        fig.add_trace(
            go.Indicator(
                mode="gauge+number+delta",
                value=metrics.ts_tv_ratio,
                domain={'x': [0, 1], 'y': [0, 1]},
                title={'text': "Ti/Tv Ratio"},
                gauge={
                    'axis': {'range': [None, 4]},
                    'bar': {'color': "darkblue"},
                    'steps': [
                        {'range': [0, 2], 'color': "lightgray"},
                        {'range': [2, 2.5], 'color': "yellow"},
                        {'range': [2.5, 4], 'color': "green"}
                    ],
                    'threshold': {
                        'line': {'color': "red", 'width': 4},
                        'thickness': 0.75,
                        'value': 2.1
                    }
                }
            ),
            row=1, col=1
        )
        
        # Ti/Tv bar chart
        fig.add_trace(
            go.Bar(x=['Transitions', 'Transversions'],
                   y=[metrics.transitions, metrics.transversions],
                   marker_color=['blue', 'red']),
            row=1, col=2
        )
        
        fig.update_layout(height=400, title_text="Ti/Tv Analysis")
        return fig
    

    @staticmethod
    def _create_position_heatmap(df: pd.DataFrame):
        """Create position density heatmap"""
        import numpy as np
        import pandas as pd
        import plotly.express as px

        # Create position bins for heatmap
        df_sample = df.sample(n=min(10000, len(df))) if len(df) > 10000 else df

        # Group by chromosome and create position bins
        heatmap_data = []
        chrom_labels = []

        for chrom in sorted(df_sample['CHROM'].unique()):
            chrom_data = df_sample[df_sample['CHROM'] == chrom]
            if len(chrom_data) > 0:
                pos_bins = pd.cut(chrom_data['POS'], bins=50, labels=False)
                bin_counts = pd.Series(pos_bins).value_counts().sort_index()
                heatmap_data.append(bin_counts.values)
                chrom_labels.append(str(chrom))

        if heatmap_data:
            # Pad each row to equal length with 0s
            max_len = max(len(row) for row in heatmap_data)
            padded_data = [
                np.pad(row, (0, max_len - len(row)), constant_values=0)
                for row in heatmap_data
            ]
            heatmap_array = np.array(padded_data)

            fig = px.imshow(
                heatmap_array,
                labels=dict(x="Position Bins", y="Chromosome", color="Variant Count"),
                x=[f"Bin {i+1}" for i in range(max_len)],
                y=chrom_labels,
                title="Variant Density Heatmap by Chromosome and Position",
                color_continuous_scale='viridis'
            )
            fig.update_layout(height=600)
            return fig

        return px.scatter(title="No data available for heatmap")

    
    @staticmethod
    def _create_quality_metrics_plot(df: pd.DataFrame):
        """Create advanced quality metrics visualization"""
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Quality by Chromosome', 'Quality by Variant Type',
                          'Quality Cumulative Distribution', 'Quality Outlier Analysis')
        )
        
        # Quality by chromosome
        for chrom in sorted(df['CHROM'].unique())[:10]:  # Limit to first 10 chromosomes
            chrom_data = df[df['CHROM'] == chrom]
            fig.add_trace(
                go.Box(y=chrom_data['QUAL'], name=chrom),
                row=1, col=1
            )
        
        # Quality by variant type
        ref_lens = df['REF'].str.len()
        alt_lens = df['ALT'].str.len()
        
        snp_mask = (ref_lens == 1) & (alt_lens == 1)
        indel_mask = (ref_lens != alt_lens)
        
        fig.add_trace(
            go.Box(y=df[snp_mask]['QUAL'], name='SNPs'),
            row=1, col=2
        )
        fig.add_trace(
            go.Box(y=df[indel_mask]['QUAL'], name='Indels'),
            row=1, col=2
        )
        
        # Cumulative distribution
        qual_sorted = np.sort(df['QUAL'])
        cumulative = np.arange(1, len(qual_sorted) + 1) / len(qual_sorted)
        
        fig.add_trace(
            go.Scatter(x=qual_sorted, y=cumulative, mode='lines',
                      name='Cumulative Distribution'),
            row=2, col=1
        )
        
        # Outlier analysis
        Q1 = df['QUAL'].quantile(0.25)
        Q3 = df['QUAL'].quantile(0.75)
        IQR = Q3 - Q1
        outliers = df[(df['QUAL'] < Q1 - 1.5 * IQR) | (df['QUAL'] > Q3 + 1.5 * IQR)]
        
        fig.add_trace(
            go.Scatter(x=outliers.index, y=outliers['QUAL'], mode='markers',
                      name='Outliers', marker=dict(color='red')),
            row=2, col=2
        )
        
        fig.update_layout(height=800, title_text="Advanced Quality Metrics Analysis")
        return fig
    
    @staticmethod
    def _create_consequence_distribution(annotated_df: pd.DataFrame):
        """Create consequence type distribution"""
        if 'consequence' not in annotated_df.columns:
            return px.bar(title="No consequence data available")
        
        consequence_counts = annotated_df['consequence'].value_counts()
        
        fig = px.bar(
            x=consequence_counts.values,
            y=consequence_counts.index,
            orientation='h',
            title="Variant Consequence Distribution",
            labels={'x': 'Number of Variants', 'y': 'Consequence Type'},
            color=consequence_counts.values,
            color_continuous_scale='plasma'
        )
        
        fig.update_layout(height=600, yaxis={'categoryorder': 'total ascending'})
        return fig
    
    @staticmethod
    def _create_impact_analysis(annotated_df: pd.DataFrame):
        """Create impact analysis visualization"""
        if 'impact' not in annotated_df.columns:
            return px.pie(title="No impact data available")
        
        impact_counts = annotated_df['impact'].value_counts()
        
        # Define impact severity colors
        impact_colors = {
            'HIGH': '#d62728',      # Red
            'MODERATE': '#ff7f0e',  # Orange
            'LOW': '#2ca02c',       # Green
            'MODIFIER': '#1f77b4'   # Blue
        }
        
        colors = [impact_colors.get(impact, '#gray') for impact in impact_counts.index]
        
        fig = px.pie(
            values=impact_counts.values,
            names=impact_counts.index,
            title="Variant Impact Distribution",
            color_discrete_map=impact_colors
        )
        
        fig.update_traces(textposition='inside', textinfo='percent+label')
        return fig
    
    @staticmethod
    def _create_frequency_analysis(annotated_df: pd.DataFrame):
        """Create allele frequency analysis"""
        if 'allele_frequency' not in annotated_df.columns:
            return px.histogram(pd.DataFrame({'x': []}), x='x', title="No frequency data available")
        
        freq_data = annotated_df[annotated_df['allele_frequency'] > 0]['allele_frequency']
        
        if len(freq_data) == 0:
            return px.histogram(pd.DataFrame({'x': []}), x='x', title="No frequency data available")
        
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=('Allele Frequency Distribution', 'Allele Frequency (Log Scale)')
        )
        
        # Linear scale
        fig.add_trace(
            go.Histogram(x=freq_data, nbinsx=50, name='Allele Frequency'),
            row=1, col=1
        )
        
        # Log scale
        fig.add_trace(
            go.Histogram(x=np.log10(freq_data + 1e-6), nbinsx=50, name='Log Allele Frequency'),
            row=2, col=1
        )
        
        fig.update_layout(height=600, title_text="Allele Frequency Analysis")
        return fig
    
    @staticmethod
    def _create_pathogenicity_scores(annotated_df: pd.DataFrame):
        """Create pathogenicity scores visualization"""
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('CADD Scores', 'SIFT Scores', 'PolyPhen Scores', 'Score Correlation')
        )
        
        # CADD scores
        if 'cadd_score' in annotated_df.columns:
            cadd_data = annotated_df[annotated_df['cadd_score'] > 0]['cadd_score']
            if len(cadd_data) > 0:
                fig.add_trace(
                    go.Histogram(x=cadd_data, nbinsx=30, name='CADD'),
                    row=1, col=1
                )
        
        # SIFT scores
        if 'sift_score' in annotated_df.columns:
            sift_data = annotated_df[annotated_df['sift_score'] > 0]['sift_score']
            if len(sift_data) > 0:
                fig.add_trace(
                    go.Histogram(x=sift_data, nbinsx=20, name='SIFT'),
                    row=1, col=2
                )
        
        # PolyPhen scores
        if 'polyphen_score' in annotated_df.columns:
            polyphen_data = annotated_df[annotated_df['polyphen_score'] > 0]['polyphen_score']
            if len(polyphen_data) > 0:
                fig.add_trace(
                    go.Histogram(x=polyphen_data, nbinsx=20, name='PolyPhen'),
                    row=2, col=1
                )
        
        # Score correlation
        if all(col in annotated_df.columns for col in ['cadd_score', 'sift_score']):
            correlation_data = annotated_df[(annotated_df['cadd_score'] > 0) & 
                                          (annotated_df['sift_score'] > 0)]
            if len(correlation_data) > 0:
                fig.add_trace(
                    go.Scatter(x=correlation_data['cadd_score'], 
                             y=correlation_data['sift_score'],
                             mode='markers', name='CADD vs SIFT'),
                    row=2, col=2
                )
        
        fig.update_layout(height=800, title_text="Pathogenicity Score Analysis")
        return fig

class AdvancedExporter:
    """Export analysis results in various formats"""
    
    @staticmethod
    def export_to_excel(df: pd.DataFrame, metrics: VCFMetrics, annotated_df: pd.DataFrame = None) -> bytes:
        """Export comprehensive analysis to Excel"""
        output = io.BytesIO()
        
        with pd.ExcelWriter(output, engine='openpyxl') as writer:
            # Main data
            df.to_excel(writer, sheet_name='Variants', index=False)
            
            # Metrics summary
            metrics_data = {
                'Metric': ['Total Variants', 'Chromosomes', 'SNPs', 'Indels', 
                          'Insertions', 'Deletions', 'Transitions', 'Transversions',
                          'Ti/Tv Ratio', 'Mean Quality', 'Median Quality', 'Quality Std'],
                'Value': [metrics.total_variants, metrics.chromosomes, metrics.snps,
                         metrics.indels, metrics.insertions, metrics.deletions,
                         metrics.transitions, metrics.transversions, metrics.ts_tv_ratio,
                         metrics.qual_mean, metrics.qual_median, metrics.qual_std]
            }
            
            pd.DataFrame(metrics_data).to_excel(writer, sheet_name='Summary', index=False)
            
            # Annotated data if available
            if annotated_df is not None and not annotated_df.empty:
                annotated_df.to_excel(writer, sheet_name='Annotated_Variants', index=False)
        
        output.seek(0)
        return output.getvalue()
    
    @staticmethod
    def export_to_csv(df: pd.DataFrame) -> str:
        """Export variants to CSV format"""
        return df.to_csv(index=False)
    
    @staticmethod
    def generate_report(df: pd.DataFrame, metrics: VCFMetrics, 
                       annotated_df: pd.DataFrame = None) -> str:
        """Generate comprehensive analysis report"""
        report = f"""
# VCF Analysis Report
Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## Summary Statistics
- **Total Variants**: {metrics.total_variants:,}
- **Chromosomes**: {metrics.chromosomes}
- **SNPs**: {metrics.snps:,} ({metrics.snps/metrics.total_variants*100:.1f}%)
- **Indels**: {metrics.indels:,} ({metrics.indels/metrics.total_variants*100:.1f}%)
  - Insertions: {metrics.insertions:,}
  - Deletions: {metrics.deletions:,}

## Quality Metrics
- **Mean Quality Score**: {metrics.qual_mean:.2f}
- **Median Quality Score**: {metrics.qual_median:.2f}
- **Quality Standard Deviation**: {metrics.qual_std:.2f}

## Transition/Transversion Analysis
- **Transitions**: {metrics.transitions:,}
- **Transversions**: {metrics.transversions:,}
- **Ti/Tv Ratio**: {metrics.ts_tv_ratio:.2f}

## Chromosome Distribution
"""
        
        # Add chromosome distribution
        chrom_counts = df['CHROM'].value_counts().sort_index()
        for chrom, count in chrom_counts.items():
            report += f"- **{chrom}**: {count:,} variants\n"
        
        if annotated_df is not None and not annotated_df.empty:
            report += f"""
## Annotation Summary
- **Annotated Variants**: {len(annotated_df):,}
- **Annotation Source**: {annotated_df['annotation_source'].iloc[0] if 'annotation_source' in annotated_df.columns else 'Unknown'}

### Consequence Distribution
"""
            if 'consequence' in annotated_df.columns:
                consequence_counts = annotated_df['consequence'].value_counts()
                for consequence, count in consequence_counts.items():
                    report += f"- **{consequence}**: {count:,} variants\n"
            
            if 'impact' in annotated_df.columns:
                report += "\n### Impact Distribution\n"
                impact_counts = annotated_df['impact'].value_counts()
                for impact, count in impact_counts.items():
                    report += f"- **{impact}**: {count:,} variants\n"
        
        return report

def main():
    """Main Streamlit application"""
    
    # Initialize session state
    if 'vcf_data' not in st.session_state:
        st.session_state.vcf_data = None
    if 'analysis_complete' not in st.session_state:
        st.session_state.analysis_complete = False
    if 'annotated_data' not in st.session_state:
        st.session_state.annotated_data = None
    
    # Header
    st.title("🧬 Advanced VCF Analysis Suite v3.0")
    st.markdown("### Real Annotation Edition")
    st.markdown("---")
    
    # Sidebar configuration
    with st.sidebar:
        st.header("⚙️ Configuration")
        
        # File upload options
        st.subheader("📁 Data Input")
        use_sample = st.checkbox("Use Sample Data", value=False)
        
        max_variants = st.slider(
            "Max Variants to Process",
            min_value=100,
            max_value=10000,
            value=1000,
            step=100,
            help="Limit processing for performance"
        )
        
        # Annotation settings
        st.subheader("🏷️ Annotation Settings")
        annotation_type = st.selectbox(
            "Annotation Method",
            ["none", "demo", "myvariant", "ensembl", "combined"],
            help="Choose annotation method (demo for testing, real APIs for production)"
        )
        
        max_annotation_variants = st.slider(
            "Max Variants to Annotate",
            min_value=10,
            max_value=200,
            value=50,
            step=10,
            help="Limit annotations due to API rate limits"
        )
        
        # Analysis options
        st.subheader("📊 Analysis Options")
        show_advanced_plots = st.checkbox("Show Advanced Plots", value=True)
        enable_export = st.checkbox("Enable Export Options", value=True)
    
    # Main content area
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.header("📤 VCF File Upload")
        
        if use_sample:
            st.info("Using sample VCF data for demonstration")
            if st.button("Generate Sample Data"):
                sample_generator = SampleDataGenerator()
                sample_vcf = sample_generator.generate_sample_vcf()
                
                with st.spinner("Processing sample VCF..."):
                    parser = OptimizedVCFParser()
                    header_lines, samples, df = parser.parse_vcf_fast(sample_vcf, max_variants)
                    
                    if not df.empty:
                        st.session_state.vcf_data = {
                            'header': header_lines,
                            'samples': samples,
                            'dataframe': df
                        }
                        st.session_state.analysis_complete = True
                        st.success(f"Sample data loaded: {len(df)} variants")
                    else:
                        st.error("Failed to generate sample data")
        else:
            uploaded_file = st.file_uploader(
                "Choose a VCF file",
                type=['vcf', 'vcf.gz'],
                help="Upload your VCF file for analysis"
            )
            
            if uploaded_file is not None:
                with st.spinner("Processing VCF file..."):
                    try:
                        # Handle gzipped files
                        if uploaded_file.name.endswith('.gz'):
                            content = gzip.decompress(uploaded_file.read()).decode('utf-8')
                        else:
                            content = uploaded_file.read().decode('utf-8')
                        
                        parser = OptimizedVCFParser()
                        header_lines, samples, df = parser.parse_vcf_fast(content, max_variants)
                        if not df.empty:
                            st.session_state.vcf_data = {
                                'header': header_lines,
                                'samples': samples,
                                'dataframe': df
                            }
                            st.session_state.analysis_complete = True
                            st.success(f"VCF loaded successfully: {len(df)} variants")
                        else:
                            st.error("No variants found in VCF file")
                    
                    except Exception as e:
                        st.error(f"Error processing VCF file: {str(e)}")
                        logger.error(f"VCF processing error: {e}")
    
    with col2:
        st.header("📋 File Info")
        if st.session_state.vcf_data:
            info_data = st.session_state.vcf_data
            st.info(f"""
            **Variants**: {len(info_data['dataframe'])}
            **Samples**: {len(info_data['samples'])}
            **Chromosomes**: {info_data['dataframe']['CHROM'].nunique()}
            """)
    
    # Analysis Section
    if st.session_state.analysis_complete and st.session_state.vcf_data:
        st.markdown("---")
        st.header("📊 Analysis Results")
        
        df = st.session_state.vcf_data['dataframe']
        
        # Generate metrics
        with st.spinner("Calculating metrics..."):
            metrics = FastStatsCalculator.calculate_comprehensive_stats(df)

        
        # Display key metrics
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Variants", f"{metrics.total_variants:,}")
        with col2:
            st.metric("SNPs", f"{metrics.snps:,}")
        with col3:
            st.metric("Indels", f"{metrics.indels:,}")
        with col4:
            st.metric("Ti/Tv Ratio", f"{metrics.ts_tv_ratio:.2f}")
        
        # Annotation Section
        if annotation_type != "none":
            st.markdown("---")
            st.header("🏷️ Variant Annotation")
            
            if annotation_type == "demo":
                st.info("Using demo annotation (simulated results)")
            else:
                st.warning(f"Using {annotation_type.upper()} annotation - API rate limits apply")
            
            if st.button("Start Annotation", key="annotate_btn"):
                with st.spinner(f"Annotating variants using {annotation_type}..."):
                    try:
                        annotator = AdvancedAnnotationManager()
                        
                        if annotation_type == "demo":
                            annotated_df = annotator.annotate_variants(
                                df, "demo", max_annotation_variants
                            )
                        elif annotation_type == "myvariant":
                            annotated_df = annotator.annotate_variants(
                                df, "myvariant", max_annotation_variants
                            )
                        elif annotation_type == "ensembl":
                            annotated_df = annotator.annotate_variants(
                                df, "ensembl", max_annotation_variants
                            )
                        elif annotation_type == "combined":
                            annotated_df = annotator.annotate_variants(
                                df, "combined", max_annotation_variants
                            )
                        
                        st.session_state.annotated_data = annotated_df
                        st.success(f"Annotation complete: {len(annotated_df)} variants annotated")
                        
                    except Exception as e:
                        st.error(f"Annotation failed: {str(e)}")
                        logger.error(f"Annotation error: {e}")
        
        # Display annotated data
        if st.session_state.annotated_data is not None:
            st.subheader("📋 Annotated Variants")
            
            # Show annotation summary
            annotated_df = st.session_state.annotated_data
            
            if 'consequence' in annotated_df.columns:
                consequence_counts = annotated_df['consequence'].value_counts()
                st.write("**Consequence Types:**")
                for consequence, count in consequence_counts.head(10).items():
                    st.write(f"- {consequence}: {count}")
            
            # Display data table
            st.dataframe(
                annotated_df.head(100),
                use_container_width=True,
                height=400
            )
        
            # Advanced Plots Section
            if show_advanced_plots:
                st.markdown("---")
                st.header("📈 Advanced Visualizations")

                with st.spinner("Generating advanced plots..."):
                    plots = AdvancedPlotGenerator.create_comprehensive_dashboard(
                        df, metrics, st.session_state.annotated_data
                    )

            
            # Display plots in tabs
            plot_tabs = st.tabs([
                "Chromosome Distribution", "Quality Analysis", "Variant Types", 
                "Ti/Tv Analysis", "Position Heatmap", "Quality Metrics"
            ])
            
            with plot_tabs[0]:
                if 'chrom_dist' in plots:
                    st.plotly_chart(plots['chrom_dist'], use_container_width=True)
            
            with plot_tabs[1]:
                if 'qual_dist' in plots:
                    st.plotly_chart(plots['qual_dist'], use_container_width=True)
            
            with plot_tabs[2]:
                if 'type_dist' in plots:
                    st.plotly_chart(plots['type_dist'], use_container_width=True)
            
            with plot_tabs[3]:
                if 'titv' in plots:
                    st.plotly_chart(plots['titv'], use_container_width=True)
            
            with plot_tabs[4]:
                if 'pos_heatmap' in plots:
                    st.plotly_chart(plots['pos_heatmap'], use_container_width=True)
            
            with plot_tabs[5]:
                if 'qual_metrics' in plots:
                    st.plotly_chart(plots['qual_metrics'], use_container_width=True)
            
            # Annotation-specific plots
            if st.session_state.annotated_data is not None:
                annotation_tabs = st.tabs([
                    "Consequences", "Impact Analysis", "Frequency", "Pathogenicity"
                ])
                
                with annotation_tabs[0]:
                    if 'consequence_dist' in plots:
                        st.plotly_chart(plots['consequence_dist'], use_container_width=True)
                
                with annotation_tabs[1]:
                    if 'impact_analysis' in plots:
                        st.plotly_chart(plots['impact_analysis'], use_container_width=True)
                
                with annotation_tabs[2]:
                    if 'frequency_analysis' in plots:
                        st.plotly_chart(plots['frequency_analysis'], use_container_width=True)
                
                with annotation_tabs[3]:
                    if 'pathogenicity_scores' in plots:
                        st.plotly_chart(plots['pathogenicity_scores'], use_container_width=True)
        
        # Export Section
        if enable_export:
            st.markdown("---")
            st.header("💾 Export Results")
            
            export_col1, export_col2, export_col3 = st.columns(3)
            
            with export_col1:
                if st.button("📊 Export to Excel"):
                    with st.spinner("Generating Excel file..."):
                        exporter = AdvancedExporter()
                        excel_data = exporter.export_to_excel(
                            df, metrics, st.session_state.annotated_data
                        )
                        
                        st.download_button(
                            label="Download Excel File",
                            data=excel_data,
                            file_name=f"vcf_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                        )
            
            with export_col2:
                if st.button("📄 Export to CSV"):
                    exporter = AdvancedExporter()
                    csv_data = exporter.export_to_csv(df)
                    
                    st.download_button(
                        label="Download CSV File",
                        data=csv_data,
                        file_name=f"vcf_variants_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                        mime="text/csv"
                    )
            
            with export_col3:
                if st.button("📋 Generate Report"):
                    with st.spinner("Generating comprehensive report..."):
                        exporter = AdvancedExporter()
                        report = exporter.generate_report(
                            df, metrics, st.session_state.annotated_data
                        )
                        
                        st.download_button(
                            label="Download Report",
                            data=report,
                            file_name=f"vcf_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.md",
                            mime="text/markdown"
                        )
        
        # Raw Data Section
        with st.expander("🔍 View Raw Data", expanded=False):
            st.subheader("VCF Header")
            header_text = "\n".join(st.session_state.vcf_data['header'][:20])
            st.text(header_text + ("..." if len(st.session_state.vcf_data['header']) > 20 else ""))
            
            st.subheader("Variant Data")
            st.dataframe(df.head(1000), use_container_width=True)
            
            st.subheader("Data Info")
            buffer = io.StringIO()
            df.info(buf=buffer)
            st.text(buffer.getvalue())
    
    # Footer
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: #666;'>
        <p>🧬 Advanced VCF Analysis Suite v3.0 | Built with Streamlit & Real Genomics APIs</p>
        <p>⚡ Features: Real-time annotation, Interactive plots, Comprehensive metrics</p>
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
