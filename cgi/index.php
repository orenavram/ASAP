<?
    include ("templates/definitions.tpl");
    include ("templates/header.tpl");
    include ("templates/help.tpl"); //defines all the help boxes' text

?>
<body>
	<form name="ASAP_form"
	      action="/cgi-bin/test.py"
	      method="post"
	      ENCTYPE="multipart/form-data"
	      onsubmit="return validate_input()"
	>
	<input type="hidden" name="example_page" value="no"> <!--------
 <font color=red>
 <center><h2>*** IMPORTANT NOTICE ***</h2></center></font>
 <center><h3>May 7th 2017: Please note, due to technical issues with our cluster system, the submission of computations to the queue might take longer than usual.<br>
 <br />We apologize for the inconvenience.<br>
 </h3></center>
 </font>
 --------->
	<p>
		<b>Sequence data</b> <font size=-1>(<A href="https://en.wikipedia.org/wiki/FASTQ_format">FASTQ format</A> only)</font>
		<ul><b>Run #1 data</b>
			<UL>
			R1 file: <input type="file" name="Run1_seq_File_R1"> <font size=-1>(<a href="example/run1/242_R1.fastq">example file</a>)</font><br>
			R2 file: <input type="file" name="Run1_seq_File_R2"> <font size=-1>(<a href="example/run1/242_R2.fastq">example file</a>)</font><br>
			</UL>
		</ul>
		<ul>
			<div id="add_run2">
				<a href="javascript:show_div_name('run2');hide_div_name('add_run2');">Add 2nd run</a><br>
			</div>
		</ul>	
		<div id="run2">
			<ul><b>Run #2 data</b>:
				<UL>
				R1 file: <input type="file" name="Run2_seq_File_R1"> <font size=-1>(<a href="example/run2/242_R1.fastq">example file</a>)</font><br>
				R2 file: <input type="file" name="Run2_seq_File_R2"> <font size=-1>(<a href="example/run2/242_R2.fastq">example file</a>)</font><br>
				</UL><br>
				<div id="add_run3">
					<a href="javascript:show_div_name('run3');hide_div_name('add_run3');">Add 3rd run</a><br>
				</div>
			</ul>
		</div>
		<div id="run3">
			<ul><b>Run #3 data</b>:
				<UL>
				R1 file: <input type="file" name="Run3_seq_File_R1"> <font size=-1>(<a href="example/run2/242_R1.fastq">example file</a>)</font><br>
				R2 file: <input type="file" name="Run3_seq_File_R2"> <font size=-1>(<a href="example/run2/242_R2.fastq">example file</a>)</font><br>
				</UL>
			</ul>
		</div>	
	</p><br><br>
	<p>
		<b>Chains: </b><input type="checkbox" name="chains" value="IGH" checked>IGH</input> <input type="checkbox" name="chains" value="IGK">IGK</input> <input type="checkbox" name="chains" value="IGL">IGL</input>
	</p>
	<p>
		<b>Organism: </b>
		<select name="MMU">
				<option value="False" selected>Human
				<option value="True">Mouse
		</select>
	</p>
<!--
	<p>
		<b>Cell counts: </b><textarea name="CellCounts" rows=1 cols=5 id="CellCounts"></textarea>
	</p>	
-->	
	<br><br>
	<div class="mC">
		<span class="mH" onclick="toggleMenu('Expert_users')">(+) Advanced options</span>
		<div id='Expert_users' class="mL"><spanid='Expert_users'>
			<p>
				MUT_BY_CDR3:&nbsp;&nbsp;<textarea name="MUT_BY_CDR3" rows=1 cols=5 id="MUT_BY_CDR3">True</textarea>
				&nbsp;&nbsp;<a href="javascript:toggle_help('help_MUT_BY_CDR3');" title="click for help"><img src=i_quest.jpg border="0"></a>
				<div id="help_MUT_BY_CDR3">
					<p>
						TBD
					</p>
				</div>
			</p>
			<p>
				Minimal threshold of reads' length (nucleotides):&nbsp;&nbsp;<textarea name="len_threshold" rows=1 cols=5 id="len_threshold">300</textarea>
				&nbsp;&nbsp;<a href="javascript:toggle_help('help_len_threshold');" title="click for help"><img src=i_quest.jpg border="0"></a>
				<div id="help_len_threshold">
					<p>
						TBD
					</p>
				</div>
			</p>
			<p>
				Minimal threshold of reads' average quality:&nbsp;&nbsp;<textarea name="qlty_threshold" rows=1 cols=5 id="qlty_threshold">20</textarea>
				&nbsp;&nbsp;<a href="javascript:toggle_help('help_qlty_threshold');" title="click for help"><img src=i_quest.jpg border="0"></a>
				<div id="help_qlty_threshold">
					<p>
						TBD
					</p>
				</div>
			</p>
			<p>
				Maximal threshold for CDR3 counts analysis in the overlapping runs:&nbsp;&nbsp;<textarea name="max_threshold" rows=1 cols=5 id="max_threshold">8</textarea>
				&nbsp;&nbsp;<a href="javascript:toggle_help('help_max_threshold');" title="click for help"><img src=i_quest.jpg border="0"></a>
				<div id="help_max_threshold">
					<p>
						TBD
					</p>
				</div>
			</p>
<!--
			<p>
				Weight score of sequences that their CDR3 count is greater than threshold:&nbsp;&nbsp;<textarea name="ws_seq" rows=1 cols=5 id="ws_seq">0.8</textarea>
				&nbsp;&nbsp;<a href="javascript:toggle_help('help_ws_seq');" title="click for help"><img src=i_quest.jpg border="0"></a>
				<div id="help_ws_seq">
					<p>
						TBD
					</p>
				</div>
			</p>
			<p>
				Weight score of CDR3s that their count is greater than threshold:&nbsp;&nbsp;<textarea name="ws_cdr" rows=1 cols=5 id="ws_cdr">0.3</textarea>
				&nbsp;&nbsp;<a href="javascript:toggle_help('help_ws_cdr');" title="click for help"><img src=i_quest.jpg border="0"></a>
				<div id="help_ws_cdr">
					<p>
						TBD
					</p>
				</div>
			</p>
			<p>
				Weight score of overlapping CDR3:&nbsp;&nbsp;<textarea name="ws_OL_cdr" rows=1 cols=5 id="ws_OL_cdr">1</textarea>
				&nbsp;&nbsp;<a href="javascript:toggle_help('help_ws_OL_cdr');" title="click for help"><img src=i_quest.jpg border="0"></a>
				<div id="help_ws_OL_cdr">
					<p>
						TBD
					</p>
				</div>
			</p>
			<p>
				Weight score of CDR3 out of cell count summary:&nbsp;&nbsp;<textarea name="ws_cc" rows=1 cols=5 id="ws_cc">0.1</textarea>
				&nbsp;&nbsp;<a href="javascript:toggle_help('help_ws_cc');" title="click for help"><img src=i_quest.jpg border="0"></a>
				<div id="help_ws_cc">
					<p>
						TBD
					</p>
				</div>
			</p>
-->			
	</div>  <!-- END  Expert_users div -->
	
	<p><br><br></p>

	<p align=left><input type="checkbox" name="send_user" value="yes" onclick="javascript:if (ASAP_form.send_user.checked) {show_div_name('user_add');}else {hide_div_name('user_add');}">Send a link to the results by e-mail (optional)
	<div id="user_add">
		<p align=left>Please enter your email address&nbsp;&nbsp;&nbsp;&nbsp;</p>
		<INPUT TYPE="TEXT" NAME="email_add" SIZE=50 bgcolor="yellow"><br>
		<font size=2>Your email address will be used to update you the moment the results are ready.</font>
	</div>
			
	<div id="confirm" style="display:none">Please confirm your email address<INPUT TYPE="HIDDEN" NAME="confirm_email_add" SIZE=50 bgcolor="yellow" ></div>
	<P align=left><br>Job title (optional)
			<INPUT TYPE="TEXT" NAME="Job_Title_txt" SIZE=50 bgcolor="yellow"><br>
			<font size=2>Enter a descriptive job title for your ASAP query</font>
	</P>		
	<br><br><br>
<!--	
	<div class="g-recaptcha" data-sitekey="6LfNOC8UAAAAAO2gQKB2YGqhURSiwZqU1avn5hHb"></div>
		<noscript>
			<div style="width: 302px; height: 352px;">
			<div style="width: 302px; height: 352px; position: relative;">
			<div style="width: 302px; height: 352px; position: absolute;">
            <iframe src="https://www.google.com/recaptcha/api/fallback?k=6LfNOC8UAAAAAO2gQKB2YGqhURSiwZqU1avn5hHb"
                    frameborder="0" scrolling="no"
                    style="width: 302px; height:352px; border-style: none;">
            </iframe>
            </div>
            <div style="width: 250px; height: 80px; position: absolute; border-style: none;
                        bottom: 21px; left: 25px; margin: 0px; padding: 0px; right: 25px;">
            <textarea id="g-recaptcha-response" name="g-recaptcha-response"
                      class="g-recaptcha-response"
                      style="width: 250px; height: 80px; border: 1px solid #c1c1c1;
                      margin: 0px; padding: 0px; resize: none;" value="">
            </textarea>
            </div>
            </div>
            </div>
        </noscript>
-->
	<br>	
	<input type=submit value="Submit" style="color:red;font-weight:bold;height:35px;font-size:15px">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
	<input type="reset" value="Clear" id="buttons" style="font-weight:bold;height:35px;font-size:15px" onclick="javascript:hide_div_name('run2');show_div_name('add_run2');javascript:hide_div_name('run3');javascript:show_div_name('add_run3');">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
	<button type="button" onclick="javascript:LoadExample()" style="font-weight:bold;height:35px;font-size:15px">Load Example</button>
	
	<br><br><br><br><br>


	
	</form>
</body>
</html>
