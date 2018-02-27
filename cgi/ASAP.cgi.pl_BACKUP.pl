#!/usr/bin/perl -w
use CGI;
# use CGI::Carp qw(fatalsToBrowser); # for debug - print fatal errors to the web-site screen
use Scalar::Util qw(looks_like_number);
use LWP::Simple;
use lib "/bioseq/bioSequence_scripts_and_constants";
use lib "/bioseq/asap";
use GENERAL_CONSTANTS;
use BIOSEQUENCE_FUNCTIONS;
use ASAP_CONSTANTS;
use Storable;
use File::Copy;
use File::Path qw(make_path);
use Excel::Writer::XLSX;

###### this command = take the name of the path of files from the HTML
my $query = new CGI;

my %VARS=();
my %FORM=();

my $stored_data_file = "input.data"; # a stored hash of %VARS
my $stored_form_data = "form.data"; # a stored hash of %FORM

my $qsub_script = "qsub.sh"; # script to run in the queue

my $submission_time;
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$submission_time = $hour . ':' . $min . ':' . $sec;
my $curr_time = $submission_time." $mday-".($mon+1)."-".($year+1900);

$FORM{sent_from_example}=$query->param('example_page'); # yes | no

if ($FORM{sent_from_example} eq "no") { # not submitted from an example page
	# 1st run data
	$FORM{Run1_userSeq_File_R1} = $query->param('Run1_seq_File_R1');
	$FORM{Run1_userSeq_File_R2} = $query->param('Run1_seq_File_R2');

	#2nd run data
	$FORM{Run2_userSeq_File_R1} = $query->param('Run2_seq_File_R1');
	$FORM{Run2_userSeq_File_R2} = $query->param('Run2_seq_File_R2');
}
else # sent from example page
{
	$FORM{Run1_userSeq_File_R1}="242_R1.fastq";
	$FORM{Run1_userSeq_File_R2}="242_R2.fastq";
	$FORM{Run2_userSeq_File_R1}="242_R1.fastq";
	$FORM{Run2_userSeq_File_R2}="242_R2.fastq";
}

$FORM{MMU} = $query->param('MMU'); # False=Human; True=mouse
$FORM{chains} = [$query->param('chains')]; # {IGH,IGK,IGL}

$FORM{raw_data_file_suffix}= "txt";

# $FORM{CellCounts}=$query->param('CellCounts');
# $FORM{CellCounts}=9999;



# advanced
$FORM{MUT_BY_CDR3} = $query->param('MUT_BY_CDR3'); # True/False
$FORM{len_threshold} = $query->param('len_threshold');
$FORM{qlty_threshold} = $query->param('qlty_threshold');
$FORM{max_threshold} = $query->param('max_threshold');
$FORM{ws_seq} = $query->param('ws_seq');
$FORM{ws_cdr} = $query->param('ws_cdr');
$FORM{ws_OL_cdr} = $query->param('ws_OL_cdr');
$FORM{ws_cc} = $query->param('ws_cc');

$FORM{userJob_Title} = $query->param('Job_Title_txt');

$FORM{user_email} = $query->param('email_add');
$FORM{email_checkbox}= $query->param('send_user');
$FORM{SPAMMERS_BOTS} = $query->param('confirm_email_add'); # This is hidden field that only spammers bots might fill in, if it is contain a value it is a spammer...


if (!defined $FORM{SPAMMERS_BOTS}) {$FORM{SPAMMERS_BOTS}="yes";}
if (!defined $FORM{user_email}){$FORM{user_email}="";}
if (!defined $FORM{email_checkbox}){$FORM{email_checkbox}="no";}
if (!defined $FORM{userJob_Title}) {$FORM{userJob_Title}="";}

# default values - NOT REALY USED ANYMORE BUT NEEDED FOR THE PARAM FILE
if (!defined $FORM{ws_seq}) {$FORM{ws_seq}=0.8;}
if (!defined $FORM{ws_cdr}) {$FORM{ws_cdr}=0.3;}
if (!defined $FORM{ws_OL_cdr}) {$FORM{ws_OL_cdr}=1;}
if (!defined $FORM{ws_cc}) {$FORM{ws_cc}=0.1;}


# ##### GOOGLE RECAPTCHA VALIDATION START######
 # my $captcha_response=$query->param('g-recaptcha-response');
 # my $secret_captcha="6LfNOC8UAAAAAMgePorCoWiWZuXtzZ3_ya0Aj8lc";
 # if ((!defined $captcha_response) or ($captcha_response eq ""))
 # {
	 # print "Content-type: text/html\n\n";
	## Note there is a newline between
	## this header and Data
   # print "<html> <head>\n";
	 # print "<title>ASAP server</title>";
	 # print "</head>\n";
	 # print "<body>\n";
	 # print "<font size=+3 color='red'>ERROR! ASAP session has been terminated: </font><br>\n";
	 # print "The re-captcha box 'I'm not a robot' was not checked, please check and submit your run again or <a href=\"mailto:bioSequence\@tauex.tau.ac.il?subject=ASAP ReCaptcha problem\">contact us</a> for assitance<br>\n";
	 # print "</body> </html>\n";
	 # die;
 # }
 # else
 # {
	 # my $URL="https://www.google.com/recaptcha/api/siteverify?secret=$secret_captcha&response=$captcha_response";
	 # my $contents = get $URL;
	 # if ($contents!~/\"success\": true/)
	 # {
		 # print "Content-type: text/html\n\n";

		## Note there is a newline between
		## this header and Data
		 # print "<html> <head>\n";
		 # print "<title>ASAP server</title>";
		 # print "</head>\n";
		 # print "<body>\n";
		 # print "<font size=+3 color='red'>ERROR!  ASAP session has been terminated: </font><br>\n";
		 # print "The re-captcha challenge has failed; Please submit your run again or <a href=\"mailto:bioSequence\@tauex.tau.ac.il?subject=ASAP ReCaptcha problem\">contact us</a> for assitance<br>\n";
		 # print "</body> </html>\n";
		 # die;
	 # }
 # }

# ##### GOOGLE RECAPTCHA VALIDATION END######

### Sending e-mails vars
$VARS{send_email_dir} = GENERAL_CONSTANTS::SEND_EMAIL_DIR_IBIS;
$VARS{smtp_server} = GENERAL_CONSTANTS::SMTP_SERVER;
$VARS{userName} = GENERAL_CONSTANTS::ADMIN_USER_NAME;
$VARS{userPass} = GENERAL_CONSTANTS::ADMIN_PASSWORD;
$VARS{IsSPAM}=0;
if ($FORM{SPAMMERS_BOTS} eq "yes"){$VARS{IsSPAM}=1;}
###### This part creates a folder for each executation (same as static in c++)
$VARS{run_number} = $^T;
my $random_number = int(rand(9999));
$VARS{run_number}.=$random_number;

###### resultsDir is where the results will be.
my $resultsDir = ASAP_CONSTANTS::ASAP_RESULTS_DIR; #XMXMXMXMX
my $logs_dir = ASAP_CONSTANTS::ASAP_LOGS_DIR;
$VARS{OutLogFile} = $logs_dir.$VARS{run_number}.".log";

#start LOG
&start_log();
open (my $LOG, ">>","$VARS{OutLogFile}");
print $LOG "The value of FORM{SPAMMERS_BOTS} is $FORM{SPAMMERS_BOTS}\n";

###### WorkingDir is where the results of the current run will be.
$VARS{WorkingDir} = $resultsDir . $VARS{run_number} . "/";
my $ans=system ("mkdir $VARS{WorkingDir}");
print $LOG "system (mkdir $VARS{WorkingDir}):$ans\n";
$ans=system ("chmod oug+wrx $VARS{WorkingDir}");
print $LOG "system (chmod oug+wrx $VARS{WorkingDir}):$ans\n";
$VARS{output_page} = "output.php";
chmod 0755, $VARS{WorkingDir};
chmod 0755, "$VARS{WorkingDir}$VARS{output_page}";

my $sample_name="reads"; # constant sample name for the webserver. ignored in the webserver. should be number otherwise
$VARS{Run1_reads_dir}=$VARS{WorkingDir}."$sample_name/run1/";
make_path($VARS{Run1_reads_dir});
if (($FORM{Run2_userSeq_File_R1} ne "") and ($FORM{Run2_userSeq_File_R2} ne ""))
{
	$VARS{Run2_reads_dir}=$VARS{WorkingDir}."$sample_name/run2/";
	make_path($VARS{Run2_reads_dir});
}

$VARS{OutputDir}=$VARS{WorkingDir}."outputs/";
# save email
if ($FORM{user_email} ne "")
{
	open (my $USER_MAIL, ">",$VARS{WorkingDir}."user_email.txt");
	print $USER_MAIL $FORM{user_email};
	close ($USER_MAIL);
	chmod 0600, $VARS{WorkingDir}."user_email.txt";
}
# create cell counts file
#$VARS{cell_counts_xlsx_file_name}="h_Illumina_summary.xlsx";
#print $LOG "make_cell_reads_xslx($VARS{WorkingDir}$VARS{cell_counts_xlsx_file_name},$sample_name,$FORM{CellCounts});\n";
#make_cell_reads_xslx("$VARS{WorkingDir}$VARS{cell_counts_xlsx_file_name}",$sample_name,$FORM{CellCounts});

# save job title
if ($FORM{userJob_Title} ne "")
{
	open (my $USER_JOB_TITLE, ">",$VARS{WorkingDir}."job_title.txt");
	print $USER_JOB_TITLE $FORM{userJob_Title};
	close ($USER_JOB_TITLE);
	chmod 0600, $VARS{WorkingDir}."job_title.txt";
}

# INPUT FILES
$VARS{fastq_R1} = "R1.fastq"; #generic fixed name for the seq file for all run #. the user file will be copied to this file. NOTE: sample name must be number!!
$VARS{fastq_R2} = "R2.fastq";



# determine and validate files suffix
$VARS{R1_R1_suffix}="";
$VARS{R1_R2_suffix}="";
$VARS{R2_R1_suffix}="";
$VARS{R2_R2_suffix}="";

if ($FORM{Run1_userSeq_File_R1}=~/(\.[^.]+)$/)
{
 $VARS{R1_R1_suffix}=$1 if ($1 ne ".fastq");
}
if ($FORM{Run1_userSeq_File_R2}=~/(\.[^.]+)$/)
{
 $VARS{R1_R2_suffix}=$1 if ($1 ne ".fastq");
}
if ($FORM{Run2_userSeq_File_R1}=~/(\.[^.]+)$/)
{
 $VARS{R2_R1_suffix}=$1 if ($1 ne ".fastq");
}
if ($FORM{Run2_userSeq_File_R2}=~/(\.[^.]+)$/)
{
 $VARS{R2_R2_suffix}=$1 if ($1 ne ".fastq");
}
# ParamFile
$VARS{ParamFile} = "parameters.txt"; # The ParamFile for the Run

###### this is the name of the program wrapper.
my $run_calc_script="/bioseq/asap/NGS_analyzer/NGS_analyzer.py";

###### here we set the html output file (where links to all files will be)
my $results_url=ASAP_CONSTANTS::ASAP_RESULTS_URL;
$VARS{run_url}=$results_url.$VARS{run_number}."/";
$VARS{output_page} = "output.php";
my $server_url=GENERAL_CONSTANTS::ASAP_URL;
###### here we set the reload interval (in seconds).

my $reload_interval = 30;

###### here we set the email of the server - for problems...
my $DEVELOPER_MAIL = GENERAL_CONSTANTS::ADMIN_EMAIL;
my $mail = "\"mailto:$DEVELOPER_MAIL?subject=ASAP%20Run%20No.:%20$VARS{run_number}\"";

###### here we set the error definitions.

my $ErrorDef = "<font size=+3 color='red'>ERROR! ASAP session has been terminated: </font>";
my $SysErrorDef = "<p><font size=+3 color='red'>SYSTEM ERROR - ASAP session has been terminated!</font><br><b>Please wait for a while and try to run ASAP again</b></p>\n";
my $ContactDef = "\n<H3><center>For assistance please <a href=$mail>contact us</a> and mention this number: $VARS{run_number}</H3>\n";
&start_output_html();

###### Move directly to the output file
print "Location: $VARS{run_url}$VARS{output_page}\n\n";
close(STDIN);
close(STDOUT);
close(STDERR);

if (($VARS{R1_R1_suffix} ne "") and ($VARS{R1_R1_suffix} ne ".gz"))
{
	exit_on_error("user_error", "Run #1 R1 file is seems to be in '$VARS{R1_R1_suffix}' format which is currently not suppurted. Please make sure it is either fastq (unzipped) file or fastq in a <a href='https://en.wikipedia.org/wiki/Gzip'>gzip</a> (.gz) format\n",$VARS{IsSPAM})
}
if (($VARS{R1_R2_suffix} ne "") and ($VARS{R1_R2_suffix} ne ".gz"))
{
	exit_on_error("user_error", "Run #1 R2 file is seems to be in '$VARS{R1_R2_suffix}' format which is currently not suppurted. Please make sure it is either fastq (unzipped) file or fastq in a <a href='https://en.wikipedia.org/wiki/Gzip'>gzip</a> (.gz) format\n",$VARS{IsSPAM})
}
if (($VARS{R2_R1_suffix} ne "") and ($VARS{R2_R1_suffix} ne ".gz"))
{
	exit_on_error("user_error", "Run #2 R1 file is seems to be in '$VARS{R2_R1_suffix}' format which is currently not suppurted. Please make sure it is either fastq (unzipped) file or fastq in a <a href='https://en.wikipedia.org/wiki/Gzip'>gzip</a> (.gz) format\n",$VARS{IsSPAM})
}
if (($VARS{R2_R2_suffix} ne "") and ($VARS{R2_R2_suffix} ne ".gz"))
{
	exit_on_error("user_error", "Run #2 R2 file is seems to be in '$VARS{R2_R2_suffix}' format which is currently not suppurted. Please make sure it is either fastq (unzipped) file or fastq in a <a href='https://en.wikipedia.org/wiki/Gzip'>gzip</a> (.gz) format\n",$VARS{IsSPAM})
}
&upload_user_data();

# TO DO: validate inputs
#if (!looks_like_number($FORM{CellCounts})) {exit_on_error ('user_error',"Cell counts must be a number...\n",$VARS{IsSPAM});}
#if ($FORM{CellCounts}<1){exit_on_error('user_error',"Cell counts must be greater than 0...\n",$VARS{IsSPAM});}

if ($FORM{sent_from_example} eq "no")
{
	if ($FORM{Run1_userSeq_File_R1} eq "") {exit_on_error('user_error',"R1 file for 'Run #1' was not provided\n",$VARS{IsSPAM});}
	if ($FORM{Run1_userSeq_File_R2} eq "") {exit_on_error('user_error',"R2 file for 'Run #1' was not provided\n",$VARS{IsSPAM});}
#	if ($FORM{Run2_userSeq_File_R1} eq "") {exit_on_error('user_error',"R1 file for 'Run #2' was not provided\n",$VARS{IsSPAM});}
#	if ($FORM{Run2_userSeq_File_R2} eq "") {exit_on_error('user_error',"R2 file for 'Run #2' was not provided\n",$VARS{IsSPAM});}
}

my $errors = "";
if ($errors ne "")
{
	exit_on_error("user_error",$errors, $VARS{IsSPAM});
}

&prepare_ParamFile(); 
&store_data();
&Update_Users_Log();
&submit_job_to_Q();
&report_new_job_submitted();
close $LOG;






# ======================= SUBRUTINES =======================
sub prepare_ParamFile
{
	open (my $PARAMS,">","$VARS{WorkingDir}$VARS{ParamFile}") || exit_on_error('sys_error',"Can't open PARAMS '$VARS{WorkingDir}$VARS{ParamFile}' $!");
	my $mixer_path=ASAP_CONSTANTS::MiXCR_dir;
	my $chains_str=join(",",@{$FORM{chains}});
	my $fastq_path="$VARS{Run1_reads_dir}";
	my $out_path="$VARS{OutputDir}run1";
	if (($FORM{Run2_userSeq_File_R1} ne "") and ($FORM{Run2_userSeq_File_R2} ne ""))
	{
		$fastq_path.=",$VARS{Run2_reads_dir}";
		$out_path.=",$VARS{OutputDir}run2";
	}
	print $PARAMS <<EndOfPARAMS;
#String that represents the working directory (run_i should be there under ruin folder; the output will be there under 'outputs' folder)
#/bioseq/asap/
$VARS{WorkingDir}

#String that represents the path to MiXCR executable file
$mixer_path

#String that represents a sample number. Should match the FASTQ file name
$sample_name

#Integer that represents number of runs (the output will be +1 because of the joint)
2

#List of strings that represent the chains
#comma delimited (if more than one) without spaces!
$chains_str

#Integer that represents minimal threshold of reads' length (nucleotides)
$FORM{len_threshold}

#Integer that represents minimal threshold of reads' average quality
$FORM{qlty_threshold}

#Boolean that indicates whether the samples originated in mice (and not human) i.e. True=Mouse; False=Human
$FORM{qlty_threshold}

#String that represent the raw data files suffix. txt / xls / etc...
$FORM{raw_data_file_suffix}

EndOfPARAMS
	close ($PARAMS);
}

sub start_log
{
    ###### open the log file
    system 'echo "(touch '.$VARS{OutLogFile}.'; chmod oug+w '.$VARS{OutLogFile}.')" | /bin/tcsh';
    open (my $LOG, ">",$VARS{OutLogFile});
    print $LOG "\n************** LOG FILE *****************\n\n";
    print $LOG "Begin time: ".BIOSEQUENCE_FUNCTIONS::printTime()."\n";
    if (($FORM{user_email} ne "") and ($FORM{email_checkbox} eq "yes"))
	{
        print $LOG "User email is: $FORM{user_email}\n";
    }
    elsif (($FORM{user_email} ne "") and ($FORM{email_checkbox} ne "yes")) 
	{
		print $LOG "user gave email but the FORM{email_checkbox} is '$FORM{email_checkbox}'\n";
    }
    else 
	{
        print $LOG "User did not provide an email address\n";
    }
	print $LOG "Run Parameters:\n";
    close $LOG;
}

###### Start writing the output web page of ASAP
sub start_output_html 
{
    my $OutHtmlFile="$VARS{WorkingDir}$VARS{output_page}";
	
	print $LOG "start_output_html: Opening the file $OutHtmlFile, and change the permissions of the WorkingDir\n";
	system 'echo "(touch '."$VARS{WorkingDir}$VARS{output_page}".'; chmod 0755 '.$OutHtmlFile.')" | /bin/tcsh';
	open (my $OUTPUT, ">","$OutHtmlFile") || exit_on_error('sys_error', "start_output_html: Cannot open the output HTML file '$OutHtmlFile' for writing",$VARS{IsSPAM});
	my $html_dir_path=ASAP_CONSTANTS::ASAP_html_dir;
	print $OUTPUT <<EndOfHTML;
<?php
\$path = \"$html_dir_path\";
set_include_path(get_include_path() . PATH_SEPARATOR . \$path);
include ("templates/definitions.tpl"); 
?>
<HTML>
<HEAD>
<META HTTP-EQUIV="REFRESH" CONTENT=$reload_interval> </HEAD>
<META HTTP-EQUIV="PRAGMA" CONTENT="NO-CACHE"> </HEAD>
<meta http-equiv=\"X-UA-Compatible\" content=\"IE=EmulateIE7\"/>
<TITLE>ASAP Run - $FORM{usrSeq_File} $FORM{userJob_Title} Number: $VARS{run_number} </TITLE>
<link rel="icon" href="/ASAP_icon.gif">

<style type="text/css">
#menu {

text-decoration: none;
        color: white;
font-size: 12px;
font-weight: 700;
}

ul.in_progress {
    list-style-image: url('$server_url/inprogress.gif');
    padding-bottom: 0 em;
}

ul.finished {
    list-style-image: url('$server_url/finished.gif');
}
</style>
<link rel="stylesheet" type="text/css" href="$server_url/ASAP.css">
<script src="$server_url/clmenu.js" type="text/javascript"></script> 
<link href="$server_url/clmenu.css" type="text/css" rel="stylesheet" /> 

</HEAD>
<H1 align=center>ASAP Job Status - <FONT color='red'>RUNNING</FONT></h1>

<blockquote>

<p><font face=Verdana>
ASAP is now processing your request.<br>

This page will be automatically updated every 30 seconds. You can also reload it manually.<br>
Once the job has finished, several links to the output files will appear below.
<br><br>
If you wish to view these results at a later time without recalculating
them, please bookmark this page. The results will be kept in the server for three months.

</font></p>

<font face=Verdana><u><h4>Running Parameters:</h4></u></font>
	
EndOfHTML

	print $OUTPUT "<font face=Verdana>\n";
	if ($FORM{userJob_Title} ne "")
	{
		print $OUTPUT "Job title: <b>$FORM{userJob_Title}</b><br><br>\n";
	}
	
	print $OUTPUT "<b>Run 1 sequence data:</b><br>\n";
	print $OUTPUT "<ul>R1 = <A HREF='$sample_name/run1/".$VARS{fastq_R1}."' TARGET=_blank>".$FORM{Run1_userSeq_File_R1}."</A><br>\n";
	print $OUTPUT "R2 = <A HREF='$sample_name/run1/".$VARS{fastq_R2}."' TARGET=_blank>".$FORM{Run1_userSeq_File_R2}."</A><br></ul>\n";
	if (($FORM{Run2_userSeq_File_R1} ne "") and ($FORM{Run2_userSeq_File_R2} ne ""))
	{
		print $OUTPUT "<b>Run 2 sequence data:</b><br>\n";
		print $OUTPUT "<ul>R1 = <A HREF='$sample_name/run2/".$VARS{fastq_R1}."' TARGET=_blank>".$FORM{Run2_userSeq_File_R1}."</A><br>\n";
		print $OUTPUT "R2 = <A HREF='$sample_name/run2/".$VARS{fastq_R2}."' TARGET=_blank>".$FORM{Run2_userSeq_File_R2}."</A><br></ul>\n";
	}

	print $OUTPUT "Organism: ";
	if ($FORM{MMU} eq "False"){print $OUTPUT "Human<br>\n";}
	else {print $OUTPUT "Mouse<br>\n";}
	print $OUTPUT "Chains: ",join(",",@{$FORM{chains}}),"<br>\n";
	# print $OUTPUT "Cell counts: ",$FORM{CellCounts},"<br>\n";
	
	
	print $OUTPUT "</font><br>\n";
	
	close ($OUTPUT);
}
sub upload_user_data{
    if ($FORM{sent_from_example} eq "yes")
	{
		# copy example data
		my $Run1_R1_file=ASAP_CONSTANTS::ASAP_html_dir.ASAP_CONSTANTS::example_file_Run1_R1;
		my $Run1_R2_file=ASAP_CONSTANTS::ASAP_html_dir.ASAP_CONSTANTS::example_file_Run1_R2;
		my $Run2_R1_file=ASAP_CONSTANTS::ASAP_html_dir.ASAP_CONSTANTS::example_file_Run2_R1;
		my $Run2_R2_file=ASAP_CONSTANTS::ASAP_html_dir.ASAP_CONSTANTS::example_file_Run2_R2;
		print $LOG "upload_file : copy example file '$Run1_R1_file' and save it as '$VARS{Run1_reads_dir}$VARS{fastq_R1}'\n";
		print $LOG "upload_file : copy example file '$Run1_R2_file' and save it as '$VARS{Run1_reads_dir}$VARS{fastq_R2}'\n";
		print $LOG "upload_file : copy example file '$Run2_R1_file' and save it as '$VARS{Run2_reads_dir}$VARS{fastq_R1}'\n";
		print $LOG "upload_file : copy example file '$Run2_R2_file' and save it as '$VARS{Run2_reads_dir}$VARS{fastq_R2}'\n";
		copy ($Run1_R1_file,"$VARS{Run1_reads_dir}$VARS{fastq_R1}");
		copy ($Run1_R2_file,"$VARS{Run1_reads_dir}$VARS{fastq_R2}");
		copy ($Run2_R1_file,"$VARS{Run2_reads_dir}$VARS{fastq_R1}");
		copy ($Run2_R2_file,"$VARS{Run2_reads_dir}$VARS{fastq_R2}");
	}
	else
	{
		# upload Seqs file

		if ((!-e "$VARS{Run1_reads_dir}$VARS{fastq_R1}$VARS{R1_R1_suffix}") and ($FORM{Run1_userSeq_File_R1} ne "")){
			upload_file($FORM{Run1_userSeq_File_R1},"$VARS{Run1_reads_dir}$VARS{fastq_R1}$VARS{R1_R1_suffix}");
		}
		if ((!-e "$VARS{Run1_reads_dir}$VARS{fastq_R2}$VARS{R1_R2_suffix}") and ($FORM{Run1_userSeq_File_R2} ne "")){
			upload_file($FORM{Run1_userSeq_File_R2},"$VARS{Run1_reads_dir}$VARS{fastq_R2}$VARS{R1_R2_suffix}");
		}
		if ((!-e "$VARS{Run2_reads_dir}$VARS{fastq_R1}$VARS{R2_R1_suffix}") and ($FORM{Run2_userSeq_File_R1} ne "")){
			upload_file($FORM{Run2_userSeq_File_R1},"$VARS{Run2_reads_dir}$VARS{fastq_R1}$VARS{R2_R1_suffix}");
		}
		if ((!-e "$VARS{Run2_reads_dir}$VARS{fastq_R2}$VARS{R2_R2_suffix}") and ($FORM{Run2_userSeq_File_R2} ne "")){
			upload_file($FORM{Run2_userSeq_File_R2},"$VARS{Run2_reads_dir}$VARS{fastq_R2}$VARS{R2_R2_suffix}");
		}
	}
	
	#try to deal with spammers TO DO
	# if (-e "$VARS{WorkingDir}$VARS{userSeq_File}")
	# {
		# $VARS{IsSPAM}=CheckForSpam("$VARS{WorkingDir}$VARS{userSeq_File}",$VARS{WorkingDir});
		# if ($FORM{SPAMMERS_BOTS} eq "yes"){$VARS{IsSPAM}=1;}
	# }
	# if (-e "$VARS{WorkingDir}$VARS{userTree_File}")
	# {
		# $VARS{IsSPAM}=CheckForSpam("$VARS{WorkingDir}$VARS{userTree_File}",$VARS{WorkingDir});
		# if ($FORM{SPAMMERS_BOTS} eq "yes"){$VARS{IsSPAM}=1;}
	# }
}

#---------------------------------------------
# upload the file to the server and convert user file from dos to unix
sub upload_file{
	my $user_file = shift;
	my $file_name_on_server = shift;  # only the name! without the path

	$CGI::POST_MAX = 1024 * 1000 * ASAP_CONSTANTS::MAX_UPLOAD_FILE_SIZE;
	if (!$user_file){
		exit_on_error('user_error', "There was a problem uploading your file '$user_file'. Please note that the maximum file size is ".ASAP_CONSTANTS::MAX_UPLOAD_FILE_SIZE."G.",$VARS{IsSPAM});
	}
#	open LOG, ">>$VARS{OutLogFile}";
	print $LOG "upload_file : upload the file '$user_file' and save it as $file_name_on_server\n";
	open (UPLOADFILE, ">$file_name_on_server") or exit_on_error('sys_error', "Cannot open the file '$file_name_on_server'  for writing $!",$VARS{IsSPAM});
	binmode UPLOADFILE;
	while (<$user_file>){
		print UPLOADFILE;
    }
    close UPLOADFILE;
    if (!-e "$file_name_on_server" or -z "$file_name_on_server"){
        exit_on_error('user_error', "Cannot upload the file \"$user_file\", Please verify that the file exists and contains data.",$VARS{IsSPAM});
    }
#	close LOG;
}

#---------------------------------------------
sub exit_on_error
{
	my $which_error = shift;
	my $error_msg = shift;
	my $IsSpam=shift; # if suspected as spam don't send email to user 0 for no 1 for yes;
	my $error_definition = "<br><br><font size=+2 color='red'>ERROR! ASAP session has been terminated:</font><br />\n";
	my $syserror = "<br><br><font size=+1 color='red'>A SYSTEM ERROR OCCURRED!</font><br />Plesae try to run ASAP again in a few minutes.<br />We apologize forthe inconvenience.<br />\n";

	if ($which_error eq 'user_error'){
		# open (my $LOG, ">>","$VARS{OutLogFile}");
		print $LOG "\n\t EXIT on error:\n$error_msg\n";
		if (-e "$VARS{WorkingDir}$VARS{output_page}") # output already created
		{
			open (my $OUTPUT,">>","$VARS{WorkingDir}/$VARS{output_page}");
			print $OUTPUT $error_definition."$error_msg";
			close ($OUTPUT);
		}
		else # OUTPUT WAS NOT CREATED print $error_msg to the screen
		{
			print "Content-type: text/html\n\n";
			print "<html>\n";
			print "<head>\n";
			print "<title>ERROR has occurred</title>\n";
			print "</head>\n";
			print "<body>\n";
			print $error_definition."$error_msg";
		}
		# close $LOG;
	}
    elsif($which_error eq 'sys_error'){
		send_administrator_mail_on_error ($error_msg) if (!$IsSpam);
		# open (my $LOG, ">>","$VARS{OutLogFile}");
		print $LOG "\n$error_msg\n";
		if (-e "$VARS{WorkingDir}$VARS{output_page}") # output already created
		{
			open (my $OUTPUT,">>","$VARS{WorkingDir}/$VARS{output_page}");
			print $OUTPUT $syserror;
			close ($OUTPUT);
		}
		else # Output was not created
		{
			print "Content-type: text/html\n\n";
			print "<html>\n";
			print "<head>\n";
			print "<title>ERROR has occurred</title>\n";
			print "</head>\n";
			print "<body>\n";
			print $syserror;
		}
		#print $error_msg to the log file
		# close $LOG;
	}
	if (($FORM{user_email} ne "") and (!$IsSpam) and ($FORM{email_checkbox} eq "yes"))
	{
		send_mail_on_error();
	}
	update_output_that_run_failed();
	# open (my $LOG, ">>","$VARS{OutLogFile}");
	print $LOG "\nSUSPECTED SPAM!!! Emails were not sent!!!\n" if ($IsSpam);
	print $LOG "\nExit Time: ".(BIOSEQUENCE_FUNCTIONS::printTime)."\n";
	close $LOG;
	chmod 0755, $VARS{working_dir};
	exit;
}

########################################################################################
sub send_mail_on_error
{
	my $email_subject;
	my $HttpPath=$VARS{run_url}."/".$VARS{output_page};
	$email_subject = "'Your ASAP run $VARS{run_number} FAILED'";
	my $email_message = "'Hello,\\n\\nUnfortunately your ASAP run (number ".$VARS{run_number}.") has failed.\\nPlease have a look at ".$HttpPath." for further details\\n\\nSorry for the inconvenience\\nASAP Team'";

	my $msg = "ssh bioseq\@lecs2 \"cd $VARS{send_email_dir}; ".'./sendEmail.pl -f \'TAU BioSequence <bioSequence@tauex.tau.ac.il>\' -t \''.$FORM{user_email}.'\' -u '.$email_subject.' -xu '.$VARS{userName}.' -xp '.$VARS{userPass}.' -s '.$VARS{smtp_server}.' -m '.$email_message.'"';
	#if ($attach ne ''){$msg.=" -a $attach"; print LOG "sending $msg\n";}
	open (my $LOG, ">>","$VARS{OutLogFile}");
	print $LOG "MESSAGE:$email_message\nCOMMAND:$msg\n";
	chdir $VARS{send_email_dir};
	my $email_system_return = `$msg`;
	unless ($email_system_return =~ /successfully/)    {
		print $LOG "send_mail: The message was not sent successfully. system returned: $email_system_return\n";
	}
	close $LOG;
}

####################################################################################
sub send_administrator_mail_on_error
{
	my $message=shift;
	my $email_subject;
	$email_subject = "'System ERROR has occurred on ASAP: $VARS{run_url}'";
	my $email_message = "'Hello,\\n\\nUnfortunately a system System ERROR has occurred on ASAP: $VARS{run_url}.\\nERROR: $message.'";
	my $Admin_Email=GENERAL_CONSTANTS::ADMIN_EMAIL;
	my $msg = "ssh bioseq\@lecs2 \"cd $VARS{send_email_dir}; ".'./sendEmail.pl -f \'bioSequence@tauex.tau.ac.il\' -t \''."bioSequence\@tauex.tau.ac.il".'\' -u '.$email_subject.' -xu '.$VARS{userName}.' -xp '.$VARS{userPass}.' -s '.$VARS{smtp_server}.' -m '.$email_message.'"';
	#if ($attach ne ''){$msg.=" -a $attach"; print LOG "sending $msg\n";}
	#  print LOG "MESSAGE:$email_message\nCOMMAND:$msg\n";
	chdir $VARS{send_email_dir};
	my $email_system_return = `$msg`;
}
####################################################################################
sub CheckForSpam # check the input for the possible of spammer...
{
	my $inFile=shift;
	my $Working_Dir=shift;
	open (my $INFILE,"<",$inFile) || exit_on_error ("sys_error","CheckForSpam can't open file: '$inFile' $!",$VARS{IsSPAM});
	my $IsSpam=0;
	while ((defined (my $line=<$INFILE>)) and ($IsSpam==0))
	{
		if ($line=~/href\s*\=\s*/)
		{
			$IsSpam=1;
			last;
		}
		elsif ($line=~/http.*?\:\/\//)
		{
			$IsSpam=1;
			last;
		}
	}
	if ($IsSpam==1)
	{
		open (FLAG,">$Working_Dir"."SPAM");
		close (FLAG);
	}
	return ($IsSpam);
}
####################################################################################
sub update_output_that_run_failed
{
#	close $OUTPUT;
	# finish the output page
	open (my $IN, "<","$VARS{WorkingDir}/$VARS{output_page}");
	my @output = <$IN>;
	close ($IN);
	
	# remove the refresh commands from the output page
	open (my $OUTPUT, ">","$VARS{WorkingDir}/$VARS{output_page}");
	foreach my $line (@output)
	{
		if (($line=~/TTP-EQUIV="REFRESH"/) or ($line=~/CONTENT="NO-CACHE"/))
		{
			next;
		}
		elsif ($line=~/(.*)RUNNING(.*)/)
		{
			print $OUTPUT $1."FAILED".$2;
		}
		else {
			print $OUTPUT $line;
		}
	}
	print $OUTPUT "<h4 class=footer align=\"center\">Questions and comments are welcome! Please <span class=\"admin_link\"><a href=\"mailto:bioSequence\@tauex.tau.ac.il\?subject=ASAP\%20Run\%20Number\%20$VARS{run_number}\">contact us</a></span></h4>";
	print $OUTPUT "</body>\n";
	print $OUTPUT "</html>\n";	
	close $OUTPUT;
}
####################################################################################
sub store_data
{
	open (my $LOG, ">>","$VARS{OutLogFile}");
	print $LOG "store_data : storing hashes to files '$stored_data_file' and '$stored_form_data'\n";
	close $LOG;
	store \%VARS, "$VARS{WorkingDir}$stored_data_file";

	if ($FORM{Run1_userSeq_File_R1} =~ /([^\/]+)$/){$FORM{Run1_userSeq_File_R1}=$1;}
	if ($FORM{Run1_userSeq_File_R2} =~ /([^\/]+)$/){$FORM{Run1_userSeq_File_R2}=$1;}
	if ($FORM{Run2_userSeq_File_R1} =~ /([^\/]+)$/){$FORM{Run2_userSeq_File_R1}=$1;}
	if ($FORM{Run2_userSeq_File_R2} =~ /([^\/]+)$/){$FORM{Run2_userSeq_File_R2}=$1;}
	store \%FORM, "$VARS{WorkingDir}$stored_form_data";
	chmod 0600, "$VARS{WorkingDir}$stored_data_file";
	chmod 0600, "$VARS{WorkingDir}$stored_form_data";
}
####################################################################################
sub Update_Users_Log
{
	my $user_ip = $ENV{'REMOTE_ADDR'};  # this is one of the variables that perl's %ENV gives us
	my $rns_log=GENERAL_CONSTANTS::ASAP_LOG;
	open LIST, ">>".$rns_log;
#	flock LIST, 2;
	print LIST $curr_time." ".$VARS{run_number}." ".$user_ip." ".$FORM{user_email}."\n";
#	flock LIST, 8;
	close LIST;

#	my $email_subject = "'New ASAP Run: $VARS{run_number}'";
#	my $email_message = "'New ASAP Run started at: $VARS{run_url}'";
#	my $Admin_Email=GENERAL_CONSTANTS::ADMIN_EMAIL;
#	my $msg = $VARS{send_email_dir}.'/sendEmail.pl -f \'TAU BioSequence <bioSequence@tauex.tau.ac.il>\' -t \''."haim.ashkenazy\@gmail.com".'\' -u '.$email_subject.' -xu '.$VARS{userName}.' -xp '.$VARS{userPass}.' -s '.$VARS{smtp_server}.' -m '.$email_message;
#	system ("ssh bioseq\@lecs2 \"cd $VARS{send_email_dir}; $msg;\"");
}

####################################################################################
sub submit_job_to_Q
{
    open (my $QSUB_SH, ">","$VARS{WorkingDir}$qsub_script") or exit_on_error('sys_error', "submit_job_to_Q : cannot open the QSUB file '$VARS{WorkingDir}$qsub_script' for writing",$VARS{IsSPAM});
### FOR LECS
####################
	print $QSUB_SH '#!/bin/tcsh',"\n";
	print $QSUB_SH '#$ -N ASAP_',"$VARS{run_number}\n";
	print $QSUB_SH '#$ -S /bin/tcsh',"\n";
	print $QSUB_SH '#$ -cwd',"\n";
	print $QSUB_SH '#$ -e ',$VARS{WorkingDir},'$JOB_NAME.$JOB_ID.ER',"\n";
	print $QSUB_SH '#$ -o ',$VARS{WorkingDir},'$JOB_NAME.$JOB_ID.OU',"\n";
	print $QSUB_SH "module load python\/anaconda_python-3.5\n";
  print $QSUB_SH 'setenv PATH "/bioseq/Programs/MAFFT_7.222/installation/bin:${PATH}"'."\n";
	print $QSUB_SH "cd $VARS{WorkingDir}\n";
    print $QSUB_SH "python3 $run_calc_script $VARS{WorkingDir}$VARS{ParamFile} \"\"\n";
	close ($QSUB_SH);
	chmod 0755, "$VARS{WorkingDir}$qsub_script";
	my $q_type="bioseq";
	my $cmd = 'ssh bioseq@lecs2 "cd '.$VARS{WorkingDir}.'; qsub '."-l $q_type $VARS{WorkingDir}$qsub_script\"";

	$VARS{qsub_job_num}="NONE";
	open (my $LOG, ">>",$VARS{OutLogFile});
	print $LOG "\nsubmit_job_to_Q :\n$cmd\n";
	my $ans = `$cmd`;
	if ($ans =~ /(\d+)/)
	{
		$VARS{qsub_job_num} = $1;
	}
	print $LOG "submit_job_to_Q : job number in the queue: $VARS{qsub_job_num}\n";
	close $LOG;
	
	open (my $QSUB_NUM, ">", "$VARS{WorkingDir}QSUB_NUN");
	print $QSUB_NUM "$VARS{qsub_job_num}\n";
	close $QSUB_NUM;

}
sub make_cell_reads_xslx {
	my $out_file=shift;
	my $sample_name=shift;
	my $cell_count=shift;
	 
	my $workbook  = Excel::Writer::XLSX->new( $out_file );
	my $worksheet = $workbook->add_worksheet();
	# first line must be empty in order to barsed correctly 
	$worksheet->write( "A2", "Sample name" );
	$worksheet->write( "A3", "$sample_name" );

	$worksheet->write( "B2", "Number of cells" );
	$worksheet->write( "B3", "$cell_count" );
	 
	$workbook->close;
}

sub report_new_job_submitted
{
	if ($FORM{userJob_Title}!~/pupkolab/i)
	{
		my $email_subject = "'New ASAP Run: $VARS{run_number}'";
		my $email_message = "'New ASAP Run started at: $VARS{run_url}'";
		my $Admin_Email=GENERAL_CONSTANTS::ADMIN_EMAIL;
		my $msg = "cd $VARS{send_email_dir}; ".$VARS{send_email_dir}.'/sendEmail.pl -f \'TAU BioSequence <bioSequence@tauex.tau.ac.il>\' -t \''."haim.ashkenazy\@gmail.com; orenavram\@gmail.com".'\' -u '.$email_subject.' -xu '.$VARS{userName}.' -xp '.$VARS{userPass}.' -s '.$VARS{smtp_server}.' -m '.$email_message;
		print $LOG "MESSAGE:$email_message\nCOMMAND:$msg\n";
		chdir $VARS{send_email_dir};
		#my $email_system_return = `$msg`;
		#unless ($email_system_return =~ /successfully/)    {
		#	print $LOG "send_mail: The message was not sent successfully. system returned: $email_system_return\n";
		#Comment out to send admin an email regarding a new run that started. Commented by Oren.
		#system ("$msg;");
	}
}