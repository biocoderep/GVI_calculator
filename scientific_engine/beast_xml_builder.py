def get_dynamic_beast_blocks(iqtree_model_string):
    """Parses an IQ-TREE model string (e.g., GTR+F+I+G4) and returns BEAST2 XML blocks."""
    base_model = "JC69"
    if "GTR" in iqtree_model_string: base_model = "GTR"
    elif "TN" in iqtree_model_string: base_model = "TN93"
    elif "HKY" in iqtree_model_string: base_model = "HKY"
    
    has_gamma = "+G" in iqtree_model_string
    has_inv = "+I" in iqtree_model_string
    
    subst_xml = ""
    prior_xml = ""
    operator_xml = ""
    log_xml = ""
    
    # Base Model
    if base_model == "JC69":
        subst_xml = '<substModel id="JC69.s:alignment" spec="JukesCantor"/>'
    elif base_model == "HKY":
        subst_xml = '''<substModel id="hky.s:alignment" spec="HKY" kappa="@kappa.s:alignment">
                            <frequencies id="estimatedFreqs.s:alignment" spec="Frequencies" frequencies="@freqParameter.s:alignment"/>
                        </substModel>'''
        prior_xml += '''<prior id="KappaPrior.s:alignment" name="distribution" x="@kappa.s:alignment">
                            <LogNormal id="LogNormalPrior1" name="distr" M="1.0" S="1.25"/>
                        </prior>'''
        operator_xml += '''<operator id="KappaScaler.s:alignment" spec="ScaleOperator" parameter="@kappa.s:alignment" scaleFactor="0.5" weight="0.1"/>
                           <operator id="FrequenciesExchanger.s:alignment" spec="DeltaExchangeOperator" delta="0.01" weight="0.1" parameter="@freqParameter.s:alignment"/>'''
        log_xml += '''<log idref="kappa.s:alignment"/>
                      <log idref="freqParameter.s:alignment"/>'''
    elif base_model == "TN93":
        subst_xml = '''<substModel id="tn93.s:alignment" spec="TN93" kappa1="@kappa1.s:alignment" kappa2="@kappa2.s:alignment">
                            <frequencies id="estimatedFreqs.s:alignment" spec="Frequencies" frequencies="@freqParameter.s:alignment"/>
                        </substModel>'''
        prior_xml += '''<prior id="Kappa1Prior.s:alignment" name="distribution" x="@kappa1.s:alignment">
                            <LogNormal id="LogNormalPrior2" name="distr" M="1.0" S="1.25"/>
                        </prior>
                        <prior id="Kappa2Prior.s:alignment" name="distribution" x="@kappa2.s:alignment">
                            <LogNormal id="LogNormalPrior3" name="distr" M="1.0" S="1.25"/>
                        </prior>'''
        operator_xml += '''<operator id="Kappa1Scaler.s:alignment" spec="ScaleOperator" parameter="@kappa1.s:alignment" scaleFactor="0.5" weight="0.1"/>
                           <operator id="Kappa2Scaler.s:alignment" spec="ScaleOperator" parameter="@kappa2.s:alignment" scaleFactor="0.5" weight="0.1"/>
                           <operator id="FrequenciesExchanger.s:alignment" spec="DeltaExchangeOperator" delta="0.01" weight="0.1" parameter="@freqParameter.s:alignment"/>'''
        log_xml += '''<log idref="kappa1.s:alignment"/>
                      <log idref="kappa2.s:alignment"/>
                      <log idref="freqParameter.s:alignment"/>'''
    elif base_model == "GTR":
        subst_xml = '''<substModel id="gtr.s:alignment" spec="GTR" rateAC="@rateAC.s:alignment" rateAG="@rateAG.s:alignment" rateAT="@rateAT.s:alignment" rateCG="@rateCG.s:alignment" rateCT="@rateCT.s:alignment" rateGT="@rateGT.s:alignment">
                            <frequencies id="estimatedFreqs.s:alignment" spec="Frequencies" frequencies="@freqParameter.s:alignment"/>
                        </substModel>'''
        prior_xml += '''<prior id="RateACPrior.s:alignment" name="distribution" x="@rateAC.s:alignment"><Gamma id="Gamma.0" name="distr" alpha="0.05" beta="10.0"/></prior>
                        <prior id="RateAGPrior.s:alignment" name="distribution" x="@rateAG.s:alignment"><Gamma id="Gamma.1" name="distr" alpha="0.05" beta="20.0"/></prior>
                        <prior id="RateATPrior.s:alignment" name="distribution" x="@rateAT.s:alignment"><Gamma id="Gamma.2" name="distr" alpha="0.05" beta="10.0"/></prior>
                        <prior id="RateCGPrior.s:alignment" name="distribution" x="@rateCG.s:alignment"><Gamma id="Gamma.3" name="distr" alpha="0.05" beta="10.0"/></prior>
                        <prior id="RateCTPrior.s:alignment" name="distribution" x="@rateCT.s:alignment"><Gamma id="Gamma.4" name="distr" alpha="0.05" beta="20.0"/></prior>
                        <prior id="RateGTPrior.s:alignment" name="distribution" x="@rateGT.s:alignment"><Gamma id="Gamma.5" name="distr" alpha="0.05" beta="10.0"/></prior>'''
        operator_xml += '''<operator id="RateACScaler.s:alignment" spec="ScaleOperator" parameter="@rateAC.s:alignment" scaleFactor="0.5" weight="0.1"/>
                           <operator id="RateAGScaler.s:alignment" spec="ScaleOperator" parameter="@rateAG.s:alignment" scaleFactor="0.5" weight="0.1"/>
                           <operator id="RateATScaler.s:alignment" spec="ScaleOperator" parameter="@rateAT.s:alignment" scaleFactor="0.5" weight="0.1"/>
                           <operator id="RateCGScaler.s:alignment" spec="ScaleOperator" parameter="@rateCG.s:alignment" scaleFactor="0.5" weight="0.1"/>
                           <operator id="RateCTScaler.s:alignment" spec="ScaleOperator" parameter="@rateCT.s:alignment" scaleFactor="0.5" weight="0.1"/>
                           <operator id="RateGTScaler.s:alignment" spec="ScaleOperator" parameter="@rateGT.s:alignment" scaleFactor="0.5" weight="0.1"/>
                           <operator id="FrequenciesExchanger.s:alignment" spec="DeltaExchangeOperator" delta="0.01" weight="0.1" parameter="@freqParameter.s:alignment"/>'''
        log_xml += '''<log idref="rateAC.s:alignment"/><log idref="rateAG.s:alignment"/><log idref="rateAT.s:alignment"/><log idref="rateCG.s:alignment"/><log idref="rateCT.s:alignment"/><log idref="rateGT.s:alignment"/>
                      <log idref="freqParameter.s:alignment"/>'''

    state_xml = ""
    if base_model != "JC69":
        state_xml += '<parameter id="freqParameter.s:alignment" spec="parameter.RealParameter" dimension="4" lower="0.0" name="stateNode" upper="1.0">0.25</parameter>\n'
    if base_model == "HKY":
        state_xml += '<parameter id="kappa.s:alignment" spec="parameter.RealParameter" lower="0.0" name="stateNode">2.0</parameter>\n'
    elif base_model == "TN93":
        state_xml += '<parameter id="kappa1.s:alignment" spec="parameter.RealParameter" lower="0.0" name="stateNode">2.0</parameter>\n'
        state_xml += '<parameter id="kappa2.s:alignment" spec="parameter.RealParameter" lower="0.0" name="stateNode">2.0</parameter>\n'
    elif base_model == "GTR":
        state_xml += '<parameter id="rateAC.s:alignment" spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>\n'
        state_xml += '<parameter id="rateAG.s:alignment" spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>\n'
        state_xml += '<parameter id="rateAT.s:alignment" spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>\n'
        state_xml += '<parameter id="rateCG.s:alignment" spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>\n'
        state_xml += '<parameter id="rateCT.s:alignment" spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>\n'
        state_xml += '<parameter id="rateGT.s:alignment" spec="parameter.RealParameter" lower="0.0" name="stateNode">1.0</parameter>\n'

    gamma_param = "false"
    inv_param = "false"
    
    if has_gamma:
        gamma_param = "true"
        state_xml += '<parameter id="gammaShape.s:alignment" spec="parameter.RealParameter" name="stateNode">1.0</parameter>\n'
        prior_xml += '<prior id="GammaShapePrior.s:alignment" name="distribution" x="@gammaShape.s:alignment"><Exponential id="Exponential.0" name="distr" mean="1.0"/></prior>\n'
        operator_xml += '<operator id="gammaShapeScaler.s:alignment" spec="ScaleOperator" parameter="@gammaShape.s:alignment" scaleFactor="0.5" weight="0.1"/>\n'
        log_xml += '<log idref="gammaShape.s:alignment"/>\n'
    
    if has_inv:
        inv_param = "true"
        state_xml += '<parameter id="proportionInvariant.s:alignment" spec="parameter.RealParameter" lower="0.0" name="stateNode" upper="1.0">0.1</parameter>\n'
        prior_xml += '<prior id="PropInvPrior.s:alignment" name="distribution" x="@proportionInvariant.s:alignment"><Uniform id="Uniform.1" name="distr" upper="1.0"/></prior>\n'
        operator_xml += '<operator id="propInvScaler.s:alignment" spec="ScaleOperator" parameter="@proportionInvariant.s:alignment" scaleFactor="0.5" weight="0.1"/>\n'
        log_xml += '<log idref="proportionInvariant.s:alignment"/>\n'
        
    sitemodel_xml = f'''<siteModel id="SiteModel.s:alignment" spec="SiteModel" 
                                  gammaCategoryCount="{4 if has_gamma else 1}" 
                                  shape="{ "@gammaShape.s:alignment" if has_gamma else "1.0" }" 
                                  proportionInvariant="{ "@proportionInvariant.s:alignment" if has_inv else "0.0" }">
                            <parameter id="mutationRate.s:alignment" spec="parameter.RealParameter" estimate="false" name="mutationRate">1.0</parameter>
                            {subst_xml}
                        </siteModel>'''

    return {
        "state": state_xml,
        "sitemodel": sitemodel_xml,
        "prior": prior_xml,
        "operator": operator_xml,
        "log": log_xml
    }
