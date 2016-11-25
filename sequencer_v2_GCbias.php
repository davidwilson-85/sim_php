<style type="text/css">
<!--
body,td,th {
	font-family: Courier New, Courier, mono;
	font-size: 14px;
}
-->
</style>

<?php 
/*
HiSeq200:        PAIRED END (500 GAUSS)    FIXED LENGTH
IonProton:       SINGLE END            VARIABLE LENGTH (200 GAUSS)


                  nbr reads x read size (pb)
read depth (X) = ----------------------------
                      genome length (pb)

				
             read depth (X) x genome size (pb)
nbr reads = ----------------------------------
                      read size (pb)
					  

TO DO:					  
------

*/
ini_set('memory_limit', '2048M'); //be very careful with this configuration directive
set_time_limit(0);
$time_start = (float) array_sum(explode(' ',microtime())); //records start time
$timestamp = date("YmdHis") ."-". sprintf('%06d', mt_rand(0, 999999)); //To give read names a unique string. This allows to combine alignments from different sequencing runs
$abs_path = dirname(__FILE__); //gets the absolute file system path of this file

//Code for examining folder with template genomes
$dir = $abs_path .'\template_genomes'; //Establishes directory in which to search
$folder_contents = scandir($dir); //Creates an array with all file names
unset($folder_contents[0]); //remove hidden file from Windows system
unset($folder_contents[1]); //remove hidden file from Windows system
$filenames = array_values($folder_contents);
$total_nbr_genomes = count($filenames);
/*
echo "<pre>";
print_r($filenames);
echo "</pre>";
*/
//GET GENOME LENGTH
$first_rec_genome = $filenames[0];
$genome = file_get_contents("template_genomes\\$first_rec_genome");
$clean_genome = strtolower(preg_replace('/\s+/', '', $genome));
$genome_length = strlen($clean_genome);
$genome_length_mb = $genome_length / 1000000;

//SET PARAMETERS
$sequencer_simulated = "p"; // i = HiSeq2000; p = Ion Proton.
$read_depth = 200;
$i_fragment_size = 500; //Only for HiSe2000
$i_fragment_sd = 100; //Only for HiSe2000
$read_length = 200;
$read_length_sd = 40; //Will only be used if IonProton selected
$base_calling_error_rate = 0.01; //I normally use 0.01. Is the rate per nucleotide, not as percentage.
$base_calling_error_rate_sd = 0.01; //I normally use 0.01
$gc_bias_strength = 100; //[0 - 100]. 0 generates a random library; 100 generates a skewed library (low representation of extreme GC content sequences). Fully explained in gc_content_bias.php.

$read_depth_per_rec_genome = $read_depth / $total_nbr_genomes;
$max_num_reads = round(($read_depth * $genome_length) / $read_length); //In each strand!!!
$max_num_reads_per_rec_genome = round(($read_depth_per_rec_genome * $genome_length) / $read_length); //In each strand!!!

//FUNCTIONS TO GENERATE RANDOM NUMBERS (FOR FRAGMENT AND READ LENGTHS) FROM A NORMAL DISTRIBUTION
function random_0_1() { 
	//Auxiliary function. Returns random number with flat distribution from 0 to 1
	return (float)rand()/(float)getrandmax(); 
}

function gauss() {	
	//Auxilary vars
	$x = random_0_1();
	$y = random_0_1();
	//Two independent variables with normal distribution N(0,1)
	$u = sqrt(-2*log($x))*cos(2*pi()*$y);
	$v = sqrt(-2*log($x))*sin(2*pi()*$y);
	//I will return only one, cause only one needed
	return $u;
}

if ($sequencer_simulated == "i") {
	function gauss_illumina($i_fragment_size, $i_fragment_sd) {
		return gauss() * $i_fragment_sd + $i_fragment_size;
	}
}

if ($sequencer_simulated == "p") {
	function gauss_proton($read_length, $read_length_sd) {
		return gauss() * $read_length_sd + $read_length;
	}
}

function gauss_errors($base_calling_error_rate, $base_calling_error_rate_sd) {	
	return gauss() * $base_calling_error_rate_sd + $base_calling_error_rate;
}


//CREATE FUNCTION TO DEFINE ARRAY WITH ALL POSSIBLE MUTANT BASES TO BE USED FOR BASECALLING ERRORS
function define_mutations_array($wt_base) {
	if($wt_base == 'A') {
		$possible_err_bases = array('T', 'C', 'G');
	} else if ($wt_base == 'T') {
		$possible_err_bases = array('A', 'C', 'G');
	} else if ($wt_base == 'C') {
		$possible_err_bases = array('A', 'T', 'G');
	} else if ($wt_base == 'G') {
		$possible_err_bases = array('A', 'T', 'C');
	} else {
		$possible_err_bases = array('A', 'T', 'C', 'G');
	}
	return $possible_err_bases;
}

//OUTPUT USEFUL INFO
if ($sequencer_simulated == "i") {
	echo 'Sequencer chosen: HiSeq2000<br>';
} else if ($sequencer_simulated == "p") {
	echo 'Sequencer chosen: IonProton<br>';
}

echo 'Genome length (Mb): '. $genome_length_mb .'<br>';
echo 'Number of (pooled) template sequences: '. $total_nbr_genomes .' ('. $total_nbr_genomes / 2 .' diploid genomes)<br>';
echo 'Chosen read_depth (X): '. $read_depth .'<br>';
echo 'Mean fragment size (pb): '. $i_fragment_size.'<br>';
echo '(Mean) read length (pb): '. $read_length .'<br><br>';
echo 'Nbr reads to create: '. $max_num_reads .' ('. $max_num_reads_per_rec_genome .' per template sequence)<br><br><br>';
echo "Creating reads...<br><br>";

//SIMULATE HISEQ2000
if ($sequencer_simulated == "i") {
	//GENERATE PAIRED-END READS
	$array1 = array('a','t','c','g');
	$array2 = array('T','A','G','C');

	$file_forward = fopen("reads\HiSeq2000_reads_forward_". $read_depth ."X.fq","w");
	$file_reverse = fopen("reads\HiSeq2000_reads_reverse_". $read_depth ."X.fq","w");

	for($genome_count = 1; $genome_count <= $total_nbr_genomes;){
		$current_rec_genome = $filenames[$genome_count - 1];
		$genome = file_get_contents("template_genomes\\$current_rec_genome");
		$clean_genome = strtolower(preg_replace('/\s+/', '', $genome));
		
		for($read_nbr = 1; $read_nbr <= round($max_num_reads_per_rec_genome / 2);) { //divided by 2 because loop generated pairs of reads
			$fragment_size = round(gauss_illumina($i_fragment_size, $i_fragment_sd));
			if ($fragment_size < $read_length) {
				$fragment_size = $read_length;
			}
			$fragment_start = mt_rand(0, $genome_length - $fragment_size);
		
			$forward_read = strtoupper(substr($clean_genome, $fragment_start, $read_length));
			
			$reverse_read = strrev(substr($clean_genome, $fragment_start + $fragment_size - $read_length, $read_length));
			$compl_reverse_read = str_replace($array1, $array2, $reverse_read);
			
			//analyze reads sequence to determine GC content
			$gc_count = (substr_count($forward_read, 'G') + substr_count($forward_read, 'C') + substr_count($compl_reverse_read, 'G') + substr_count($compl_reverse_read, 'C'));
			$at_count = (substr_count($forward_read, 'A') + substr_count($forward_read, 'T') + substr_count($compl_reverse_read, 'A') + substr_count($compl_reverse_read, 'T'));
			$gc_content = $gc_count / ($gc_count + $at_count + 0.0001); //it ranges from 0 to 1
			$distance_to_neutral_gc = Abs(0.5 - $gc_content);
			$gc_random1 = (mt_rand(0,50)) / 100; //generate a float number with 2 decimals between 0 and 0.50
			$gc_random2 = mt_rand(0,100);
			if ($gc_random1 < $distance_to_neutral_gc) {	
				if ($gc_random2 >= $gc_bias_strength) {
					$gc_result = 1;
				}
			} else {
				$gc_result = 1;
			}
			
			//consider reads only if passed GC bias filter. Process read and write to .fq file
			if (isset($gc_result)) {

				//SET ABSOLUTE NUMBER OF ERRORS IN FORWARD READ
				$total_nbr_err = round(strlen($forward_read) * gauss_errors($base_calling_error_rate, $base_calling_error_rate_sd));
				if ($total_nbr_err < 1) {
					$total_nbr_err = 0;
				}

				//FOR EACH ERROR IN FORWARD READ, RANDOMLY CHOOSE A POSITION AND A NEW BASE
				$all_err_pos = array();
				for($i = 1; $i <= $total_nbr_err;) {
					$err_pos = mt_rand(0, strlen($forward_read) -1); //PHP gives 0 to first string position!
					if (!in_array($err_pos, $all_err_pos)) { //if new random number not already in array
						$all_err_pos[] = $err_pos;
						$wt_base = substr($forward_read, $err_pos, 1);
						$possible_err_bases = define_mutations_array($wt_base); //calls the function define_mutations_array() and assigns the returned array to a variable
						$err_base = mt_rand(0, 2);
						$forward_read[$err_pos] = $possible_err_bases[$err_base];
						$i++;
					}
				}
				unset($all_err_pos);
				
				//SET ABSOLUTE NUMBER OF ERRORS IN REVERSE READ
				$total_nbr_err = round(strlen($reverse_read) * gauss_errors($base_calling_error_rate, $base_calling_error_rate_sd));
				if ($total_nbr_err < 1) {
					$total_nbr_err = 0;
				}

				//FOR EACH MUTATION IN REVERSE READ, RANDOMLY CHOOSE A POSITION AND A NEW BASE
				$all_err_pos = array();
				for($i = 1; $i <= $total_nbr_err;) {
					$err_pos = mt_rand(0, strlen($compl_reverse_read) -1); //PHP gives 0 to first string position!
					if (!in_array($err_pos, $all_err_pos)) { //if new random number not already in array
						$all_err_pos[] = $err_pos;
						$wt_base = substr($compl_reverse_read, $err_pos, 1);
						$possible_err_bases = define_mutations_array($wt_base); //calls the function define_mutations_array() and assigns the returned array to a variable
						$err_base = mt_rand(0, 2); //If $wt_base == N, $err_base = mt_rand(0, 3)!!!!!!!!!!
						$compl_reverse_read[$err_pos] = $possible_err_bases[$err_base];
						$i++;
					}
				}
				unset($all_err_pos);

				$fastq_string = preg_replace('/./', 'I', $forward_read); //I = ASCII73 (value 40 in Sanger encoding)
				
				//Determine which read will be first and which second in the pair
				$pos_in_pair = mt_rand(0,1); //randomly choose between two options
				if ($pos_in_pair == 1) {
					$swap_var = $forward_read;
					$forward_read = $compl_reverse_read;
					$compl_reverse_read = $swap_var;
				}
				//echo $forward_read.'<br>'. $compl_reverse_read.'<br><br>';
				
				fwrite($file_forward, "@exp". $timestamp ."-gnm". $genome_count ."-read". $read_nbr ."/1 i=". $fragment_start ."\r\n". $forward_read ."\r\n+\r\n". $fastq_string ."\r\n");
				fwrite($file_reverse, "@exp". $timestamp ."-gnm". $genome_count ."-read". $read_nbr ."/2 i=". $fragment_start ."\r\n". $compl_reverse_read ."\r\n+\r\n". $fastq_string ."\r\n");
				
				$read_nbr ++;
			}
			
			$forward_read = null;
			$reverse_read = null;
			$compl_reverse_read = null;
			$swap_var = null;
			$fastq_string = null;
			$total_nbr_err = null;
			$all_err_pos = null;
			$i = null;
			$err_pos = null;
			$wt_base = null;
			$possible_err_bases = null;
			$err_base = null;
			$gc_count = null;
			$at_count = null;
			$gc_content = null;
			$distance_to_neutral_gc = null;
			$gc_random1 = null;
			$gc_random2 = null;
			$gc_result = null;
		}
		
		$genome_count++;
		
		$current_rec_genome = null; //I use $var = null instead of unset($var) because unset does not necessarily empty immediately the memory used by the variable.
		$genome = null;  
		$clean_genome = null;
		$read_nbr = null;
	}
	fclose($file_forward);
	fclose($file_reverse);
}

//SIMULATE IONPROTON
if ($sequencer_simulated == "p") {	
	//GENETATE READS
	$array1 = array('a','t','c','g');
	$array2 = array('T','A','G','C');
	
	$file = fopen("reads\IonProton_reads_". $read_depth ."X.fq","w");
	
	for($genome_count = 1; $genome_count <= $total_nbr_genomes;) {
		$current_rec_genome = $filenames[$genome_count - 1];
		$genome = file_get_contents("template_genomes\\$current_rec_genome");
		$clean_genome = strtolower(preg_replace('/\s+/', '', $genome));
		
		for($read_nbr = 1; $read_nbr <= $max_num_reads_per_rec_genome;) {
			$read_length_current = round(gauss_proton($read_length, $read_length_sd));
			if ($read_length_current < 5) {
				$read_length_current = 5;
			}
			$read_start = mt_rand(0, $genome_length - $read_length_current);
			$strand = mt_rand(0,1); //randomly choose between two options (forward or reverse strands)
			
			if ($strand == 0) {
				$read = strtoupper(substr($clean_genome, $read_start, $read_length_current));
				$strand = "for";
			} else if ($strand == 1) {
				$reverse_read = strrev(substr($clean_genome, $read_start, $read_length_current)); //CHECK THAT THIS IS CORRECT!!!!!!!!!!!!!!!!!!
				$read = str_replace($array1, $array2, $reverse_read);
				$strand = "rev";
			}
			
			//analyze read sequence to determine GC content
			$gc_count = substr_count($read, 'G') + substr_count($read, 'C');
			$at_count = substr_count($read, 'A') + substr_count($read, 'T');
			$gc_content = $gc_count / ($gc_count + $at_count + 0.0001); //it ranges from 0 to 1
			$distance_to_neutral_gc = Abs(0.5 - $gc_content);
			$gc_random1 = (mt_rand(0,50)) / 100; //generate a float number with 2 decimals between 0 and 0.50
			$gc_random2 = mt_rand(0,100);
			if ($gc_random1 < $distance_to_neutral_gc) {	
				if ($gc_random2 >= $gc_bias_strength) {
					$gc_result = 1;
				}
			} else {
				$gc_result = 1;
			}
			
			//consider read only if passed GC bias filter. Process read and write to .fq file
			if (isset($gc_result)) {
				
				//SET ABSOLUTE NUMBER OF ERRORS IN READ
				$total_nbr_err = round(strlen($read) * gauss_errors($base_calling_error_rate, $base_calling_error_rate_sd));
				if ($total_nbr_err < 1) {
					$total_nbr_err = 0;
				}

				//FOR EACH ERROR, RANDOMLY CHOOSE A POSITION AND A NEW BASE
				$all_err_pos = array();
				for($i = 1; $i <= $total_nbr_err;) {
					$err_pos = mt_rand(0, strlen($read) - 1); //PHP gives 0 to first string position!
					if (!in_array($err_pos, $all_err_pos)) { //if new random number not already in array
						$all_err_pos[] = $err_pos;
						$wt_base = substr($read, $err_pos, 1);
						$possible_err_bases = define_mutations_array($wt_base); //calls the function define_mutations_array() and assigns the returned array to a variable
						$err_base = mt_rand(0, 2); //If $wt_base == N, $err_base = mt_rand(0, 3)!!!!!!!!!!
						$read[$err_pos] = $possible_err_bases[$err_base];
						$i++;
					}
				}
				unset($all_err_pos);
				
				$fastq_string = preg_replace('/./', 'I', $read); //I = ASCII73 (value 40)			  
				fwrite($file, "@exp". $timestamp ."-gnm". $genome_count ."-read". $read_nbr ." s=". $strand ." i=". $read_start ." l=". $read_length_current ."\r\n". $read ."\r\n+\r\n". $fastq_string ."\r\n");
					
				$read_nbr ++;
			}
			
			$read_length_current = null;
			$read_start = null;
			$strand = null;
			$read = null;
			$reverse_read = null;
			$fastq_string = null;
			$total_nbr_err = null;
			$all_err_pos = null;
			$i = null;
			$err_pos = null;
			$wt_base = null;
			$possible_err_bases = null;
			$err_base = null;
			$gc_count = null;
			$at_count = null;
			$gc_content = null;
			$distance_to_neutral_gc = null;
			$gc_random1 = null;
			$gc_random2 = null;
			$gc_result = null;
		}
		
		$genome_count++;
		
		$current_rec_genome = null; //I use $var = null instead of unset($var) because unset does not necessarily empty immediately the memory used by the variable.
		$genome = null;  
		$clean_genome = null;
		$read_nbr = null;
	}
}

echo 'Reads created.<br><br>';
$time_end = (float) array_sum(explode(' ',microtime())); //records end time
echo "Execution time: ". sprintf("%.4f", $time_end - $time_start) ." seconds."; 

?>