import os

NEW_CONTENT = '''import os
import subprocess
import re
import logging
from datetime import datetime
from scipy.stats import linregress
from Bio import Phylo
from config import Config
from .parsers import sanitize_for_phipack
from .beast_xml_builder import get_dynamic_beast_blocks

logger = logging.getLogger(__name__)

def chunk_fasta(fasta_path, max_seqs=60):
    """Splits an aligned FASTA file into multiple chunks of max_seqs to process sequentially."""
    from Bio import SeqIO
    import random
    
    records = list(SeqIO.parse(fasta_path, "fasta"))
    random.shuffle(records)
    
    num_chunks = (len(records) + max_seqs - 1) // max_seqs
    chunk_paths = []
    for i in range(num_chunks):
        chunk_records = records[i*max_seqs : (i+1)*max_seqs]
        for rec in chunk_records:
            if ',' in rec.id:
                rec.id = rec.id.replace(',', '_')
                rec.name = rec.id
                rec.description = rec.id
        chunk_path = f"{fasta_path}.chunk{i}.fasta"
        SeqIO.write(chunk_records, chunk_path, "fasta")
        chunk_paths.append(chunk_path)
        
    return chunk_paths

def calculate_dynamic_chain_length(fasta_path):
    from Bio import SeqIO
    records = list(SeqIO.parse(fasta_path, "fasta"))
    num_seqs = len(records)
    if num_seqs == 0: return 500000
    avg_len = sum(len(r.seq) for r in records) / num_seqs
    base_chain = 2000000
    dynamic_length = int(base_chain * (num_seqs / 30.0) * (avg_len / 1700.0))
    return max(500000, min(dynamic_length, 10000000))

def get_decimal_date(date_str):
    import datetime
    from dateutil import parser
    try:
        if len(date_str) == 4 and date_str.isdigit():
            return float(date_str)
        dt = parser.parse(date_str)
        year_start = datetime.datetime(dt.year, 1, 1)
        year_end = datetime.datetime(dt.year + 1, 1, 1)
        decimal = dt.year + ((dt - year_start).total_seconds() / (year_end - year_start).total_seconds())
        return round(decimal, 4)
    except Exception:
        import re
        match = re.search(r'\\b(19|20)\\d{2}\\b', date_str)
        if match: return float(match.group(0))
        return 2020.0

def extract_dates_from_fasta(fasta_path):
    from Bio import SeqIO
    date_mapping = {}
    for rec in SeqIO.parse(fasta_path, "fasta"):
        parts = rec.description.split('|')
        date_str = parts[-1].strip() if len(parts) > 1 else rec.id.split('_')[-1]
        date_mapping[rec.id] = get_decimal_date(date_str)
    return date_mapping

def detect_best_evolutionary_model(fasta_path):
    """Runs IQ-TREE ModelFinder Plus to get the best standard substitution model."""
    cmd = [Config.IQTREE_BIN, "-s", fasta_path, "-mset", "JC,HKY,TN,GTR", "-m", "TESTONLY", "-nt", "2", "-redo"]
    if Config.USE_WSL_FALLBACK:
        cmd = ["wsl"] + cmd
        
    logger.info("Detecting best substitution model with IQ-TREE ModelFinder...")
    try:
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=300)
    except Exception as e:
        logger.error(f"Model detection failed: {e}")
        return "GTR+G" # fallback
        
    iqtree_log = fasta_path + ".iqtree"
    best_model = "GTR+G"
    if os.path.exists(iqtree_log):
        with open(iqtree_log, 'r') as f:
            content = f.read()
            match = re.search(r'Best-fit model: ([\w\+\-]+)', content)
            if match:
                best_model = match.group(1)
    
    logger.info(f"Detected Model: {best_model}")
    return best_model

def get_relaxed_clock_blocks(num_taxa):
    state = f"""<parameter id="ucldMean.c:alignment" spec="parameter.RealParameter" name="stateNode">0.001</parameter>
                <parameter id="ucldStdev.c:alignment" spec="parameter.RealParameter" lower="0.0" name="stateNode">0.1</parameter>
                <stateNode id="rateCategories.c:alignment" spec="parameter.IntegerParameter" dimension="{max(1, 2*num_taxa-2)}" name="stateNode">1</stateNode>"""
    
    prior = """<prior id="ucldMeanPrior.c:alignment" name="distribution" x="@ucldMean.c:alignment"><Uniform id="Uniform.0" name="distr" upper="Infinity"/></prior>
               <prior id="ucldStdevPrior.c:alignment" name="distribution" x="@ucldStdev.c:alignment"><Gamma id="Gamma.0s" name="distr" alpha="0.5396" beta="0.3819"/></prior>"""
    
    branch = """<branchRateModel id="RelaxedClock.c:alignment" spec="beast.evolution.branchratemodel.UCRelaxedClockModel" rateCategories="@rateCategories.c:alignment" tree="@Tree.t:alignment" clock.rate="@ucldMean.c:alignment">
                    <LogNormal id="LogNormalDistributionModel.c:alignment" S="@ucldStdev.c:alignment" meanInRealSpace="true" name="distr">
                        <parameter id="RealParameter.1" spec="parameter.RealParameter" estimate="false" lower="0.0" name="M" upper="1.0">1.0</parameter>
                    </LogNormal>
                </branchRateModel>"""
    
    operators = """<operator id="ucldMeanScaler.c:alignment" spec="ScaleOperator" parameter="@ucldMean.c:alignment" scaleFactor="0.5" weight="1.0"/>
                   <operator id="ucldStdevScaler.c:alignment" spec="ScaleOperator" parameter="@ucldStdev.c:alignment" scaleFactor="0.5" weight="3.0"/>
                   <operator id="CategoriesRandomWalk.c:alignment" spec="IntRandomWalkOperator" parameter="@rateCategories.c:alignment" weight="10.0" windowSize="1"/>
                   <operator id="CategoriesSwapOperator.c:alignment" spec="SwapOperator" intparameter="@rateCategories.c:alignment" weight="10.0"/>
                   <operator id="CategoriesUniform.c:alignment" spec="UniformOperator" parameter="@rateCategories.c:alignment" weight="10.0"/>"""
    
    logs = """<log idref="ucldMean.c:alignment"/>
              <log idref="ucldStdev.c:alignment"/>
              <log id="rateStat.c:alignment" spec="beast.evolution.branchratemodel.RateStatistic" branchratemodel="@RelaxedClock.c:alignment" tree="@Tree.t:alignment"/>"""
              
    return state, prior, branch, operators, logs

def generate_beast_xml(fasta_path, output_xml, detected_model, chain_length=2000000):
    """Generates a coalescent BEAST2 XML for recombination/clock rate estimation using detected models."""
    from Bio import SeqIO
    
    records = list(SeqIO.parse(fasta_path, "fasta"))
    num_taxa = len(records)
    
    seqs_xml = ""
    for rec in records:
        seqs_xml += f'        <sequence id="seq_{rec.id}" spec="Sequence" taxon="{rec.id}" totalcount="4" value="{str(rec.seq)}"/>\\n'
        
    date_mapping = extract_dates_from_fasta(fasta_path)
    trait_string = ",\\n".join([f"{taxon}={date}" for taxon, date in date_mapping.items()])
    
    dyn_blocks = get_dynamic_beast_blocks(detected_model)
    clk_state, clk_prior, clk_branch, clk_op, clk_log = get_relaxed_clock_blocks(num_taxa)
    
    xml_content = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<beast beautitemplate='Standard' beautistatus='' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.math.distributions:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood:beast.evolution.branchratemodel" required="CoalRe v0.0.7" version="2.6">
    <data id="alignment" spec="beast.evolution.alignment.Alignment" name="alignment">
{seqs_xml}
    </data>

    <map name="Uniform" >beast.math.distributions.Uniform</map>
    <map name="Exponential" >beast.math.distributions.Exponential</map>
    <map name="LogNormal" >beast.math.distributions.LogNormalDistributionModel</map>
    <map name="Normal" >beast.math.distributions.Normal</map>
    <map name="Beta" >beast.math.distributions.Beta</map>
    <map name="Gamma" >beast.math.distributions.Gamma</map>
    <map name="LaplaceDistribution" >beast.math.distributions.LaplaceDistribution</map>
    <map name="prior" >beast.math.distributions.Prior</map>
    <map name="InverseGamma" >beast.math.distributions.InverseGamma</map>
    <map name="OneOnX" >beast.math.distributions.OneOnX</map>

    <run id="mcmc" spec="beast.core.MCMC" chainLength="{chain_length}">
        <state id="state" spec="beast.core.State" storeEvery="1000">
            <tree id="Tree.t:alignment" spec="beast.evolution.tree.Tree" name="stateNode">
                <trait id="dateTrait.t:alignment" spec="beast.evolution.tree.TraitSet" traitname="date-forward" value="{trait_string}">
                    <taxa id="TaxonSet.alignment" spec="beast.evolution.alignment.TaxonSet">
                        <alignment idref="alignment"/>
                    </taxa>
                </trait>
                <taxonset idref="TaxonSet.alignment"/>
            </tree>
            {dyn_blocks['state']}
            {clk_state}
            <parameter id="popSize.t:alignment" spec="parameter.RealParameter" name="stateNode">0.3</parameter>
            <parameter id="recombinationRate" spec="parameter.RealParameter" name="stateNode">0.05</parameter>
        </state>

        <distribution id="posterior" spec="util.CompoundDistribution">
            <distribution id="prior" spec="util.CompoundDistribution">
                <distribution id="CoalescentConstant.t:alignment" spec="Coalescent">
                    <populationModel id="ConstantPopulation.t:alignment" spec="ConstantPopulation" popSize="@popSize.t:alignment"/>
                    <treeIntervals id="TreeIntervals.t:alignment" spec="TreeIntervals" tree="@Tree.t:alignment"/>
                </distribution>
                {dyn_blocks['prior']}
                {clk_prior}
                <prior id="PopSizePrior.t:alignment" name="distribution" x="@popSize.t:alignment">
                    <OneOnX id="OneOnX.1" name="distr"/>
                </prior>
                <prior id="RecombPrior" name="distribution" x="@recombinationRate">
                    <Exponential name="distr" mean="1.0"/>
                </prior>
            </distribution>
            <distribution id="likelihood" spec="util.CompoundDistribution" useThreads="true">
                <distribution id="treeLikelihood.alignment" spec="ThreadedTreeLikelihood" data="@alignment" tree="@Tree.t:alignment">
                    {dyn_blocks['sitemodel']}
                    {clk_branch}
                </distribution>
            </distribution>
        </distribution>

        {dyn_blocks['operator']}
        {clk_op}
        <operator id="PopSizeScaler.t:alignment" spec="ScaleOperator" parameter="@popSize.t:alignment" weight="3.0"/>
        <operator id="RecombScaler" spec="ScaleOperator" parameter="@recombinationRate" weight="1.0"/>
        <operator id="CoalescentConstantTreeScaler.t:alignment" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:alignment" weight="3.0"/>
        <operator id="CoalescentConstantTreeRootScaler.t:alignment" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:alignment" weight="3.0"/>
        <operator id="CoalescentConstantUniformOperator.t:alignment" spec="beast.evolution.operators.Uniform" tree="@Tree.t:alignment" weight="30.0"/>
        <operator id="CoalescentConstantSubtreeSlide.t:alignment" spec="SubtreeSlide" tree="@Tree.t:alignment" weight="15.0"/>
        <operator id="CoalescentConstantNarrow.t:alignment" spec="Exchange" tree="@Tree.t:alignment" weight="15.0"/>
        <operator id="CoalescentConstantWide.t:alignment" spec="Exchange" isNarrow="false" tree="@Tree.t:alignment" weight="3.0"/>
        <operator id="CoalescentConstantWilsonBalding.t:alignment" spec="WilsonBalding" tree="@Tree.t:alignment" weight="3.0"/>

        <logger id="tracelog" spec="Logger" fileName="beast_recomb.log" logEvery="1000" model="@posterior" sanitiseHeaders="true" sort="smart">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
            <log idref="treeLikelihood.alignment"/>
            <log id="TreeHeight.t:alignment" spec="beast.evolution.tree.TreeHeightLogger" tree="@Tree.t:alignment"/>
            <log idref="popSize.t:alignment"/>
            <log idref="recombinationRate"/>
            {dyn_blocks['log']}
            {clk_log}
        </logger>

        <logger id="treelog.t:alignment" spec="Logger" fileName="run.trees" logEvery="1000" mode="tree">
            <log id="TreeWithMetaDataLogger.t:alignment" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:alignment"/>
        </logger>
    </run>
</beast>"""
    
    with open(output_xml, 'w') as f:
        f.write(xml_content)

def calculate_beast_recombination(fasta_path, detected_model="GTR+G"):
    chunk_paths = chunk_fasta(fasta_path, max_seqs=60)
    recombination_rates = []
    beast_evo_rates = []
    weights = []
    beast_bin = "/mnt/e/Inlfuenza/avian_influenza/GVI_index_calculator/tools/beast/bin/beast"
    
    for i, chunk in enumerate(chunk_paths):
        from Bio import SeqIO
        seq_count = len(list(SeqIO.parse(chunk, "fasta")))
        
        chain_length = calculate_dynamic_chain_length(chunk)
        xml_path = chunk + ".beast.xml"
        generate_beast_xml(chunk, xml_path, detected_model, chain_length)
        
        log_file = "beast_recomb.log"
        chunk_recomb = 0.05
        chunk_evo = 0.001
        
        cmd = [beast_bin, "-overwrite", xml_path]
        if Config.USE_WSL_FALLBACK: cmd = ["wsl"] + cmd
            
        try:
            logger.info(f"Starting BEAST MCMC on chunk {i+1}/{len(chunk_paths)}")
            work_dir = os.path.dirname(chunk)
            subprocess.run(cmd, cwd=work_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=1800, text=True)
            
            log_path = os.path.join(work_dir, log_file)
            if os.path.exists(log_path):
                import pandas as pd
                df = pd.read_csv(log_path, sep='\\t', comment='#')
                burn_in = int(len(df) * 0.1)
                if 'recombinationRate' in df.columns:
                    valid_samples = df['recombinationRate'].iloc[burn_in:]
                    if not valid_samples.empty: chunk_recomb = float(valid_samples.mean())
                if 'ucldMean.c:alignment' in df.columns:
                    valid_samples = df['ucldMean.c:alignment'].iloc[burn_in:]
                    if not valid_samples.empty: chunk_evo = float(valid_samples.mean())
                    
                recombination_rates.append(chunk_recomb)
                beast_evo_rates.append(chunk_evo)
                weights.append(seq_count)
        except Exception as e:
            logger.error(f"BEAST chunk {i+1} failed: {e}")
        finally:
            if os.path.exists(chunk): os.remove(chunk)
            if os.path.exists(xml_path): os.remove(xml_path)
            
    avg_recomb = float(np.average(recombination_rates, weights=weights)) if recombination_rates else 0.05
    avg_evo = float(np.average(beast_evo_rates, weights=weights)) if beast_evo_rates else 0.001
    return avg_recomb, avg_evo

def generate_bdsky_xml(fasta_path, output_xml, detected_model, chain_length=2000000):
    from Bio import SeqIO
    records = list(SeqIO.parse(fasta_path, "fasta"))
    num_taxa = len(records)
    
    seqs_xml = ""
    for rec in records:
        seqs_xml += f'        <sequence id="seq_{rec.id}" spec="Sequence" taxon="{rec.id}" totalcount="4" value="{str(rec.seq)}"/>\\n'
        
    date_mapping = extract_dates_from_fasta(fasta_path)
    trait_string = ",\\n".join([f"{taxon}={date}" for taxon, date in date_mapping.items()])
    
    dyn_blocks = get_dynamic_beast_blocks(detected_model)
    clk_state, clk_prior, clk_branch, clk_op, clk_log = get_relaxed_clock_blocks(num_taxa)
    
    xml_content = f"""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<beast beautitemplate='Standard' beautistatus='' namespace="beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood:beast.evolution.branchratemodel:bdsky.evolution.speciation:bdsky.evolution.tree:bdsky.math.distributions:bdsky.util" required="BDSKY v1.4.5" version="2.6">
    <data id="alignment" spec="beast.evolution.alignment.Alignment" name="alignment">
{seqs_xml}
    </data>
    
    <map name="Uniform" >beast.math.distributions.Uniform</map>
    <map name="Exponential" >beast.math.distributions.Exponential</map>
    <map name="LogNormal" >beast.math.distributions.LogNormalDistributionModel</map>
    <map name="Normal" >beast.math.distributions.Normal</map>
    <map name="Beta" >beast.math.distributions.Beta</map>
    <map name="Gamma" >beast.math.distributions.Gamma</map>
    <map name="prior" >beast.math.distributions.Prior</map>

    <run id="mcmc" spec="beast.core.MCMC" chainLength="{chain_length}">
        <state id="state" spec="beast.core.State" storeEvery="1000">
            <tree id="Tree.t:alignment" spec="beast.evolution.tree.Tree" name="stateNode">
                <trait id="dateTrait.t:alignment" spec="beast.evolution.tree.TraitSet" traitname="date-forward" value="{trait_string}">
                    <taxa id="TaxonSet.alignment" spec="beast.evolution.alignment.TaxonSet">
                        <alignment idref="alignment"/>
                    </taxa>
                </trait>
                <taxonset idref="TaxonSet.alignment"/>
            </tree>
            {dyn_blocks['state']}
            {clk_state}
            <parameter id="reproductiveNumber.t:alignment" spec="parameter.RealParameter" lower="0.0" name="stateNode" upper="10.0">2.0</parameter>
            <parameter id="becomeUninfectiousRate.t:alignment" spec="parameter.RealParameter" lower="0.0" name="stateNode" upper="10.0">1.0</parameter>
            <parameter id="samplingProportion.t:alignment" spec="parameter.RealParameter" lower="0.0" name="stateNode" upper="1.0">0.01</parameter>
        </state>

        <distribution id="posterior" spec="util.CompoundDistribution">
            <distribution id="prior" spec="util.CompoundDistribution">
                <distribution id="BDSKY.t:alignment" spec="bdsky.evolution.speciation.BirthDeathSkylineModel" becomeUninfectiousRate="@becomeUninfectiousRate.t:alignment" reproductiveNumber="@reproductiveNumber.t:alignment" samplingProportion="@samplingProportion.t:alignment" tree="@Tree.t:alignment"/>
                {dyn_blocks['prior']}
                {clk_prior}
                <prior id="RePrior.t:alignment" name="distribution" x="@reproductiveNumber.t:alignment">
                    <LogNormal id="LogNormal.1" name="distr" meanInRealSpace="true" M="2.0" S="1.25"/>
                </prior>
                <prior id="becomeUninfectiousRatePrior.t:alignment" name="distribution" x="@becomeUninfectiousRate.t:alignment">
                    <LogNormal id="LogNormal.2" name="distr" meanInRealSpace="true" M="36.5" S="1.25"/>
                </prior>
                <prior id="samplingProportionPrior.t:alignment" name="distribution" x="@samplingProportion.t:alignment">
                    <Beta id="Beta.1" name="distr" alpha="1.0" beta="100.0"/>
                </prior>
            </distribution>
            <distribution id="likelihood" spec="util.CompoundDistribution" useThreads="true">
                <distribution id="treeLikelihood.alignment" spec="ThreadedTreeLikelihood" data="@alignment" tree="@Tree.t:alignment">
                    {dyn_blocks['sitemodel']}
                    {clk_branch}
                </distribution>
            </distribution>
        </distribution>

        {dyn_blocks['operator']}
        {clk_op}
        <operator id="ReScaler.t:alignment" spec="ScaleOperator" parameter="@reproductiveNumber.t:alignment" weight="10.0"/>
        <operator id="becomeUninfectiousRateScaler.t:alignment" spec="ScaleOperator" parameter="@becomeUninfectiousRate.t:alignment" weight="3.0"/>
        <operator id="samplingProportionScaler.t:alignment" spec="ScaleOperator" parameter="@samplingProportion.t:alignment" weight="10.0"/>
        <operator id="BDSKYTreeScaler.t:alignment" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:alignment" weight="3.0"/>
        <operator id="BDSKYTreeRootScaler.t:alignment" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:alignment" weight="3.0"/>
        <operator id="BDSKYUniformOperator.t:alignment" spec="beast.evolution.operators.Uniform" tree="@Tree.t:alignment" weight="30.0"/>
        <operator id="BDSKYSubtreeSlide.t:alignment" spec="SubtreeSlide" tree="@Tree.t:alignment" weight="15.0"/>
        <operator id="BDSKYNarrow.t:alignment" spec="Exchange" tree="@Tree.t:alignment" weight="15.0"/>
        <operator id="BDSKYWide.t:alignment" spec="Exchange" isNarrow="false" tree="@Tree.t:alignment" weight="3.0"/>
        <operator id="BDSKYWilsonBalding.t:alignment" spec="WilsonBalding" tree="@Tree.t:alignment" weight="3.0"/>

        <logger id="tracelog" spec="Logger" fileName="bdsky.log" logEvery="1000" model="@posterior" sanitiseHeaders="true" sort="smart">
            <log idref="posterior"/>
            <log idref="likelihood"/>
            <log idref="prior"/>
            <log idref="treeLikelihood.alignment"/>
            <log id="TreeHeight.t:alignment" spec="beast.evolution.tree.TreeHeightLogger" tree="@Tree.t:alignment"/>
            <log idref="reproductiveNumber.t:alignment"/>
            <log idref="becomeUninfectiousRate.t:alignment"/>
            <log idref="samplingProportion.t:alignment"/>
            {dyn_blocks['log']}
            {clk_log}
        </logger>
    </run>
</beast>"""
    
    with open(output_xml, 'w') as f:
        f.write(xml_content)

def calculate_bdsky_re(fasta_path, detected_model="GTR+G"):
    chunk_paths = chunk_fasta(fasta_path, max_seqs=60)
    re_rates = []
    weights = []
    beast_bin = "/mnt/e/Inlfuenza/avian_influenza/GVI_index_calculator/tools/beast/bin/beast"
    
    for i, chunk in enumerate(chunk_paths):
        from Bio import SeqIO
        seq_count = len(list(SeqIO.parse(chunk, "fasta")))
        
        chain_length = calculate_dynamic_chain_length(chunk)
        xml_path = chunk + ".bdsky.xml"
        generate_bdsky_xml(chunk, xml_path, detected_model, chain_length)
        
        log_file = "bdsky.log"
        chunk_re = 1.2
        
        cmd = [beast_bin, "-overwrite", xml_path]
        if Config.USE_WSL_FALLBACK: cmd = ["wsl"] + cmd
            
        try:
            logger.info(f"Starting BDSKY MCMC on chunk {i+1}/{len(chunk_paths)}")
            work_dir = os.path.dirname(chunk)
            subprocess.run(cmd, cwd=work_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=1800, text=True)
            
            log_path = os.path.join(work_dir, log_file)
            if os.path.exists(log_path):
                import pandas as pd
                df = pd.read_csv(log_path, sep='\\t', comment='#')
                if 'reproductiveNumber' in df.columns:
                    burn_in = int(len(df) * 0.1)
                    valid_samples = df['reproductiveNumber'].iloc[burn_in:]
                    if not valid_samples.empty:
                        chunk_re = float(valid_samples.mean())
                        re_rates.append(chunk_re)
                        weights.append(seq_count)
        except Exception as e:
            pass
        finally:
            if os.path.exists(chunk): os.remove(chunk)
            if os.path.exists(xml_path): os.remove(xml_path)
            
    if re_rates:
        import numpy as np
        return float(np.average(re_rates, weights=weights))
    return 1.2

def calculate_true_evolutionary_rate(fasta_path, detected_model="GTR+G"):
    """Calculates evolutionary rate using IQ-TREE with LSD2 Molecular Clock."""
    tree_newick = ""
    
    dates_dict = extract_dates_from_fasta(fasta_path)
    dates_file = fasta_path + ".dates.txt"
    with open(dates_file, "w") as f:
        f.write(f"{len(dates_dict)}\\n")
        for k, v in dates_dict.items():
            f.write(f"{k} {v}\\n")
            
    cmd = [Config.IQTREE_BIN, "-s", fasta_path, "-m", detected_model, "--date", dates_file, "-T", "2", "-redo"]
    if Config.USE_WSL_FALLBACK:
        cmd = ["wsl"] + cmd
        
    try:
        subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=600)
    except Exception as e:
        logger.error(f"IQ-TREE execution failed: {e}")
        return 0.0, tree_newick
        
    lsd_report = fasta_path + ".lsd"
    evo_rate = 0.0
    if os.path.exists(lsd_report):
        with open(lsd_report, 'r') as f:
            content = f.read()
            match = re.search(r'RATE OF EVOLUTION:\\s+([\\d\\.eE\\-\\+]+)', content, re.IGNORECASE)
            if match: evo_rate = float(match.group(1))
    
    # We want to return the raw tree for visualization
    treefile = fasta_path + ".treefile"
    if os.path.exists(treefile):
        with open(treefile, 'r') as f:
            raw_newick = f.read().strip()
            tree_newick = re.sub(r'([A-Za-z0-9_\\.]+)\\|[^:,()]+', r'\\1', raw_newick)
            
    # Cleanup
    for ext in [".treefile", ".iqtree", ".log", ".bionj", ".ckp.gz", ".mldist", ".model.gz", ".dates.txt", ".lsd", ".date.treefile", ".timetree.nex", ".timetree.lsd"]:
        cleanup_path = fasta_path + ext
        if os.path.exists(cleanup_path): os.remove(cleanup_path)
                
    return evo_rate, tree_newick
'''

with open("E:/Inlfuenza/avian_influenza/GVI_index_calculator/scientific_engine/external_tools.py", "w") as f:
    f.write(NEW_CONTENT)
