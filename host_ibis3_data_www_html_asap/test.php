<?
    include ("templates/definitions.tpl");
    include ("templates/header.tpl");
    include ("templates/help.tpl"); //defines all the help boxes' text

?>
<body>
	<form name="ASAP_form"
	      action="/cgi-bin/ASAP_cgi.py"
	      method="post"
	      ENCTYPE="multipart/form-data"
	      onsubmit="return validate_input()"
	>
	<input type="hidden" name="example_page" value="no"> 
	<p>
		<b>Sequence data</b> <font size=-1>(<A href="https://en.wikipedia.org/wiki/FASTQ_format">FASTQ format</A> only)</font>
		<ul><b>Run #1 data</b>
			<UL>
			R1 file: <input type="file" name="run1_R1"> <font size=-1>(<a href="example/run1/242_R1.fastq">example file</a>)</font><br>
			R2 file: <input type="file" name="run1_R2"> <font size=-1>(<a href="example/run1/242_R2.fastq">example file</a>)</font><br>
			</UL>
		</ul>
		<ul><b>Run #2 data</b>
			<UL>
			R1 file: <input type="file" name="run2_R1"> <font size=-1></font><br>
			R2 file: <input type="file" name="run2_R2"> <font size=-1></font><br>
			</UL>
		</ul>
		<ul><b>Run #3 data</b>
			<UL>
			R1 file: <input type="file" name="run3_R1"> <font size=-1></font><br>
			R2 file: <input type="file" name="run3_R2"> <font size=-1></font><br>
			</UL>
		</ul>
	</p><br><br>
	<p>
		<b>Chains: </b><input type="checkbox" name="chains" value="IGH" checked>IGH</input> <input type="checkbox" name="chains" value="IGK">IGK</input> <input type="checkbox" name="chains" value="IGL">IGL</input>
	</p>
	<p>
		<b>Organism: </b>
		<select name="MMU">
				<option value="Human" selected>Human
				<option value="Mouse">Mouse
		</select>
	</p>
	<br><br>
	<div class="mC">
		<span class="mH" onclick="toggleMenu('Expert_users')">(+) Advanced options</span>
		<div id='Expert_users' class="mL"><spanid='Expert_users'>
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
                Raw data format:&nbsp;&nbsp;<textarea name="raw_data_suffix" rows=1 cols=5 id="raw_data_suffix">txt</textarea>
                &nbsp;&nbsp;<a href="javascript:toggle_help('help_raw_data_suffix');" title="click for help"><img src=i_quest.jpg border="0"></a>
                <div id="help_raw_data_suffix">
                    <p>
                        TBD
                    </p>
                </div>
            </p>
            <p>
                Add mass_spec sequence to final fasta file&nbsp;&nbsp;<input type="checkbox" name="add_mass_spec_seq" value="yes" checked>
                &nbsp;&nbsp;<a href="javascript:toggle_help('help_add_mass_spec_seq');" title="click for help"><img src=i_quest.jpg border="0"></a>
                <div id="help_add_mass_spec_seq">
                    <p>
                        TBD
                    </p>
                </div>
            </p>
		</div>  <!-- END  Expert_users div -->
	</div>

	<p><br><br></p>
	<br>
	
	<p align=left><br><input type="checkbox" name="send_user" value="yes" onclick="javascript:if (ASAP_form.send_user.checked) {show_div_name('user_add');}else {hide_div_name('user_add');}">	Send a link to the results by e-mail (optional)
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
	<br><br>
	</form>
	<a href="http://asap.tau.ac.il/test.php" class="button">
		Clear
	</a>
	 &nbsp;
	<a href="http://asap.tau.ac.il/test_example.php" class="button">
		Load example
	</a>
</body>
</html>