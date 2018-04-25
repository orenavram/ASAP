<?
    include ("templates/definitions.tpl");
    include ("templates/header.tpl");
    include ("templates/help.tpl"); //defines all the help boxes' text

?>




<h1 align=center><font  size="6" COLOR="#003399">ASAP overview</font> <br></h1>

<ul shape="square">
	<li><b><a href="#Introduction">Introduction</a></li></b>
	<li><b><a href="#good_for">What is ASAP good for?</a></li></b>
	<li><b><a href="#Input">Input</a></li></b>
	<li><b><a href="#Output">Output and methodology</a></li></b>
 <!--
	<li><b><a href="#Methodology">Methodology</a></b></li>
	<ul>
		<li><a href="#summ_stats">SpartaABC's summary statistics</a></li>
		<li><a href="#Euclidean">Euclidean distance</a></li>
		<li><a href="#search_range">Parameter search range</a></li>
	</ul>
 -->
	<li><b><a href="#References">References</a></b></li>
</ul>

<br>

<a NAME="Introduction"></a>
<h3>Introduction</h3>
<p align="justify">
Next-generation sequencing (NGS) of antibody V genes has emerged as a powerful tool in systems immunology by providing quantitative molecular information on antibody polyclonal composition [<a href="#REF_1">1-3</a>]. Reproducible and robust information on antibody repertoires is valuable for basic and applied immunology studies. However, major computational challenges exist when analyzing antibody sequences. For example, handling errors that are introduced during sample preparation and sequencing, which leads to markedly inaccurate measurements of antibody diversity.
</p>

<br><br>

<a NAME="good_for"></a>
<h3>What is ASAP good for?</h3>
<p align="justify">
ASAP server is user-friendly, free and open to all users and there is no login requirement. ASAP server generates a high confidence antibody V gene sequence archive from technical replicates generated during NGS library preparation and sequenced by NGS platforms. ASAP is applicable for researchers that are interested to address basic questions related to B cell development and differentiation as well as applied researchers that are interested in vaccine development, immunodiagnostic discovery and monoclonal antibody engineering.
</p>

<br><br>

<a NAME="Input"></a>
<h3>Input</h3>
<p align="justify">
Paired-end (R1 and R2) from one or more datasets (technical replicates) in <a href="https://en.wikipedia.org/wiki/fastq">Fastq format</a>. These datasets can be derived from high throughput sequencing of human or murine antibody V genes. Default parameters are pre-defined thus enabling the usage by non-expert users. However, advanced users can further select parameters to customize their analysis.  
</p>

<br><br>

<a NAME="Output"></a>
<h3>Output and methodology<br>TODO: ADD MORE INFO REGARDING THE OUTPUT</h3>
<p align="justify">By far the highest diversity in antibodies occurs within the CDR3 region, which is overwhelmingly responsible for antigen recognition [<a href="#REF_4">4,5</a>]. The CDR3 is thus considered as a unique identifier of antibody clones. Our methodology is divided into individual and joint parts.</p>
<p align="justify">The individual part applies standard approaches of antibody V gene sequences filtration and annotation using germline sequences from the <a href="http://www.imgt.org/">international immunogenetics information system</a> (IMGT) [<a href="#REF_6">6</a>]. The ASAP webserver utilizes MixCR [<a href="#REF_7">7</a>] to facilitate the statistical calculations of somatic hypermutation rates, CDR3 lengths, V(D)J family assignments and V(D)J recombination linkage.</p>
<p align="justify">The joint part of the ASAP webserver, comprises a novel strategy named <b>O</b>verlap <b>R</b>egion <b>D</b>etection by <b>D</b>ouble <b>S</b>equencing (ORDDS) to establish a valid reliable archive V gene sequences based on technical replicates. In ORDDS we examine the overlapping CDR3 regions and their frequency in each technical replicate to determine which sequences are less likely to be a result of a sample preparation derived and/or sequencing errors. Moreover, ORRDS enables the flexibility to define the overlap CDR3 frequency threshold in each replicate based on the statistical data generated during the run. The output of the webserver includes detailed log files, graphical representation of the analysis and a final high confidence antibody V gene sequence archive in a Fasta format. Moreover, the final Fasta file includes important metadata regarding the repertoire features for each antibody sequence that can be readily used in downstream bioinformatics analyses.	
</p>


<br><br>

<a NAME="References"></a>
<h3>References</h3>
<ol>
<li><a name="REF_1">Lavinder, J. J. et al. <i>Proc Natl Acad Sci USA</i> 111, 2259-2264 (2014).</a></li>
<br>
<li><a name="REF_2">Wine, Y. et al. <i>Proc Natl Acad Sci USA</i> 110, 2993-2998 (2013).</a></li>
<br>
<li><a name="REF_3">Georgiou, G. et al. <i>Nat Biotech</i> 32, 158-168 (2014).</a></li>
<br>
<li><a name="REF_4">Xu, J. L. et al. <i>Immunity</i> 13, 37-45 (2000).</a></li>
<br>
<li><a name="REF_5">Murphy, K. et al. <i>Janeway's immunobiology</i>.  (Garland Science, 2012).</a></li>
<br>
<li><a name="REF_6">Lefranc, M.-P. et al. <i>Nucleic Acids Research</i> 27, 209-212 (1999).</a></li>
<br>
<li><a name="REF_7">Bolotin, D. A. et al. <i>Nature Methods</i> 12, 380 (2015).</a></li>
</ol>


<BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><br><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR><BR>


</body>
</html>







