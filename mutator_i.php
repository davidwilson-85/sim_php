<?php
ini_set('memory_limit', '1024M'); //be very careful with this configuration directive
set_time_limit(0);
$abs_path = dirname(__FILE__); //gets the absolute file system path of this file

//SET PARAMETERS (i)
$total_nbr_mutations = $_GET['total_nbr_mutations'];
$mutator_mode = $_GET['mutator_mode'];

/*
//SET PARAMETERS (manually)
$total_nbr_mutations = 2;
$mutator_mode = 'li'; //Choose between d, e, li (d = genetic drift SNP mutations, e = EMS SNP mutations GC>AT, li = long insertion (t-dna, transposon)) 
*/

//Code for examining folder with template genome
$dir = $abs_path .'\input_genome'; //Establishes directory in which to search
$folder_contents = scandir($dir); //Creates an array with all file names
$file1 = $folder_contents[2];
//echo $file1 ."<br>";

//create temp file (noheader version) of template file in "temp" directory and header for output fasta file
$file1_cont = fopen("input_genome\\$file1", "r");
$temp = fopen("temp\\$file1", "w");

while(!feof($file1_cont)) {
	$line = fgets($file1_cont);
	if ($line[0] != '>') {
		fwrite($temp, $line);
	} else {
		$header = $line;
	}
}

fclose($file1_cont);
fclose($temp);

if (isset($header)) {
	$header_string = $header;
} else {
	$header_string = ">mutator_output";
}

//Code for examining temp folder with temporary template genome
$dir = $abs_path .'\temp'; //Establishes directory in which to search
$folder_contents = scandir($dir); //Creates an array with all file names
$file1_temp = $folder_contents[2];

//READ WT GENOME FILE AND REMOVE SPACES
$raw_genome = file_get_contents("temp\\$file1_temp");
$wt_genome = strtoupper(preg_replace('/\s+/', '', $raw_genome));
$genome_length = strlen($wt_genome);

//remove temp files
$temp_to_delete = glob('temp\*'); // get all file names
foreach($temp_to_delete as $delete){ // iterate files
  if(is_file($delete))
    unlink($delete); // delete file
}

//IF SNP MODES SELECTED, FOR EACH MUTATION, RANDOMLY CHOOSE A POSITION AND A NEW BASE
if ($mutator_mode == 'd' OR $mutator_mode == 'e') {
	$mut_genome = $wt_genome;
	$all_mut_pos = array();
	for($i = 1; $i <= $total_nbr_mutations;) {
		$mut_pos = mt_rand(0, $genome_length -1); //PHP gives 0 to first string position!
		if (!in_array($mut_pos, $all_mut_pos)) { //if new random number not already in array

			$wt_base = substr($mut_genome, $mut_pos, 1);
			if ($mutator_mode == 'd') { //If mode is genetic drift, take into account the four bases but not Ns
				if($wt_base == 'A') {
					$possible_mt_bases = array('T', 'C', 'G');
					$mut_base = mt_rand(0, 2);
					$all_mut_pos[] = $mut_pos;
					$mut_genome[$mut_pos] = $possible_mt_bases[$mut_base];
					$i++;
				} else if ($wt_base == 'T') {
					$possible_mt_bases = array('A', 'C', 'G');
					$mut_base = mt_rand(0, 2);
					$all_mut_pos[] = $mut_pos;
					$mut_genome[$mut_pos] = $possible_mt_bases[$mut_base];
					$i++;
				} else if ($wt_base == 'C') {
					$possible_mt_bases = array('A', 'T', 'G');
					$mut_base = mt_rand(0, 2);
					$all_mut_pos[] = $mut_pos;
					$mut_genome[$mut_pos] = $possible_mt_bases[$mut_base];
					$i++;
				} else if ($wt_base == 'G') {
					$possible_mt_bases = array('A', 'T', 'C');
					$mut_base = mt_rand(0, 2);
					$all_mut_pos[] = $mut_pos;
					$mut_genome[$mut_pos] = $possible_mt_bases[$mut_base];
					$i++;
				}
			}
			
			if ($mutator_mode == 'e') { //If mode is EMS-induced mutations, only affect G and C bases
				if ($wt_base == 'G') {
					$mut_base = 'A';
					$all_mut_pos[] = $mut_pos;
					$mut_genome[$mut_pos] = $mut_base;
					$i++;                           //Only count up valid if wt_base is G or C
				} else if ($wt_base == 'C') {
					$mut_base = 'T';
					$all_mut_pos[] = $mut_pos;
					$mut_genome[$mut_pos] = $mut_base;
					$i++;
				}
			}
		}
	}
}

//IF LARGE INSERTION FOR EACH MUTATION, RANDOMLY CHOOSE A POSITION AND A NEW BASE
if ($mutator_mode == 'li') {
	//how string positions are handled:
	// ATCG     .ATCG     A.TCG     AT.CG     ATC.G     ATCG.
	// 01234    0          1          2          3          4
	
	//Code for examining folder with insertion sequence
	$dir = $abs_path .'\input_insertion_sequence'; //Establishes directory in which to search
	$folder_contents = scandir($dir); //Creates an array with all file names
	$file2 = $folder_contents[2];
	//echo $file2 ."<br>";

	//create temp file (noheader version) of template file in "temp" directory
	$file2_cont = fopen("input_insertion_sequence\\$file2", "r");
	$temp = fopen("temp\\$file2", "w");

	while(!feof($file2_cont)) {
		$line = fgets($file2_cont);
		if ($line[0] != '>') {
			fwrite($temp, $line);
		}
	}

	fclose($file2_cont);
	fclose($temp);

	//Code for examining temp folder with insertion sequences
	$dir = $abs_path .'\temp'; //Establishes directory in which to search
	$folder_contents = scandir($dir); //Creates an array with all file names
	$file2_temp = $folder_contents[2];

	//READ WT GENOME FILE AND REMOVE SPACES
	$tdna_seq = file_get_contents("temp\\$file2_temp");
	$clean_tdna_seq = strtoupper(preg_replace('/\s+/', '', $tdna_seq));
	$tdna_length = strlen($clean_tdna_seq);
	
	$all_mut_pos = array();
	$mut_genome = '';
	for ($i = 1; $i <= $total_nbr_mutations;) {
		$mut_pos = mt_rand(1, $genome_length - 1); //To avoid insertions being added at the beginning or the end of the chromosome
		if (!in_array($mut_pos, $all_mut_pos)) { //if new random number not already in array
			//echo $mut_pos .'<br>';
			$all_mut_pos[] = $mut_pos;
			$i++;
		}
	}
	
	$ends_plus_all_mut_pos = $all_mut_pos;
	$ends_plus_all_mut_pos[] = 0;
	$ends_plus_all_mut_pos[] = $genome_length;
	sort($ends_plus_all_mut_pos);
	$ends_plus_all_mut_pos_count = count($ends_plus_all_mut_pos);

	for($fragment_count = 1; $fragment_count < $ends_plus_all_mut_pos_count;) {
		$fragment_length = $ends_plus_all_mut_pos[$fragment_count] - $ends_plus_all_mut_pos[$fragment_count - 1];
		$fragment = substr($wt_genome, $ends_plus_all_mut_pos[$fragment_count - 1], $fragment_length);
		
		//echo $fragment .' - '. $fragment_length .'<br>';
		$mut_genome .= $fragment;
		if ($fragment_count < $ends_plus_all_mut_pos_count - 1) {
			$mut_genome .= $clean_tdna_seq;
		}
		$fragment_count++;
	}
}

//CREATE A FILE WITH INFO OF EACH MUTATION CREATED
if ($mutator_mode == 'd' OR $mutator_mode == 'e') {
	sort($all_mut_pos);
	$all_mut_pos_file = fopen("output\all_mutations.txt","w");
	fwrite($all_mut_pos_file, "pos\twt_base\tmt_base\r\n");

	foreach($all_mut_pos as $pos) {
		$wt_base = substr($wt_genome, $pos, 1);
		$mut_base = substr($mut_genome, $pos, 1);
		//echo $pos + 1 ." - ". $wt_base ." - ". $mut_base ."<br>";
		fwrite($all_mut_pos_file, $pos + 1 ."\t". $wt_base ."\t". $mut_base ."\r\n");
	}
	fclose($all_mut_pos_file);
	
} else if ($mutator_mode == 'li') {
	sort($all_mut_pos);
	$all_mut_pos_file = fopen("output\all_mutations.txt","w");
	fwrite($all_mut_pos_file, "pos\r\n");

	foreach($all_mut_pos as $pos) {
		//echo $pos + 1 ."<br>";
		fwrite($all_mut_pos_file, $pos + 1 ."\r\n");
	}
	fclose($all_mut_pos_file);
}

//WRITE MUTANT SEQUENCE TO A FILE
//$mut_genome = strtolower($mut_genome); //comment this out if want to keep the mutated bases as uppercase (to easily spot them)
$mut_genome_file = fopen("output\mutated_genome.fa","w");
//fwrite($mut_genome_file, $header_string ."\r\n");
fwrite($mut_genome_file, $mut_genome);
fclose($mut_genome_file);

//remove temp files
$temp_to_delete = glob('temp\*'); // get all file names
foreach($temp_to_delete as $delete){ // iterate files
  if(is_file($delete))
    unlink($delete); // delete file
}

echo "Template: ". $file1 ."<br>Total number of mutations: ". $_GET['total_nbr_mutations'] ."<br>Mode: ". $_GET['mutator_mode'] ."<br><br>";
echo 'Done. Check "\mutator\output" folder.<br><br><a href="">Go back to index</a>';
?>
