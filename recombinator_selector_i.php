<?php
set_time_limit(10000);
ini_set('memory_limit', '512M');

//INPUT DATA
$genomeA_raw = file_get_contents('parental_genomes\ler-mut1_chr4.fa'); //if r or d mode, this file must contain mutation to map
$genomeA = preg_replace('/\s+/', '', $genomeA_raw);
$genomeB_raw = file_get_contents('parental_genomes\ler-mut2_chr4.fa');
$genomeB = preg_replace('/\s+/', '', $genomeB_raw);

//PARAMETERS
$mutation_inherintance_mode = "dr"; //"r" = recessive, "d" = dominant mt-phe, "di" = dominant wt-phe, "dr" = double recessive
$rec_freq = 0.000000061; //In chr4 simulations, values must be: $mean = 0.000000061, $sdev = 0.000000049 to fit Salome et al 2012.
$rec_freq_sd = 0.000000049;
$mutation1_position = 3068708;
$mutation2_position = 12158647; //only used in dr mode
$mutation1_genome = "genomeA"; //if "r" or "d" mode, this file must contain the mutation to be map
$mutation2_genome = "genomeB"; //only used in "dr" mode
$mutation1_color = "blue";
$mutation2_color = "red";
$total_rec_genomes = 200;

//ANALYZE INPUT DATA. CHECK THAT BOTH GENOMES ARE EQUAL IN LENGTH
$genomeA_length = strlen($genomeA);
$genomeB_length = strlen($genomeB);
if ($genomeA_length != $genomeB_length) {
	echo "Genomes lengths are not equal.<br><br>";
	exit();
} else {
	echo "genome length (Mb): ". $genomeA_length / 1000000 ."<br><br>";
}
/*
//DELETE ALL FILES FROM IMAGES FOLDER AND RECOMBINANANT GENOMES FOLDER
$files = glob('C:\wamp\www\other_files\david_wilson\NGS\recombinator\mut_rec_genomes'); //get all file names
foreach($files as $file){ //iterate files
  if(is_file($file))
    unlink($file); //delete file
}

$dir = 'C:\wamp\www\other_files\david_wilson\NGS\recombinator\mut_rec_genomes';
$dir_contents = scandir($dir);
foreach($dir_contents as $file){ //iterate files
  if(is_file($file))
    unlink($file); //delete file
}
*/

//FUNCTIONS TO CREATE GAUSSIAN DISTR. NUMBERS -> DETERMINE NBR OF CROSSOVERS FOR THE RECOMBINANT CHROMOSOMES
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
function gauss_rec_freq($rec_freq, $rec_freq_sd) {
	return gauss() * $rec_freq_sd + $rec_freq;
}

//DECLARE ARRAY TO DRAW IMAGE
$array_for_image = array();

//CREATE A NUMBER OF RECOMBINANT GENOMES FROM ISOGENIC PARENTALS A AND B
for($genome_count = 001; $genome_count <= $total_rec_genomes;){
	
	//RANDOMLY CHOOSE BETWEEN TWO OPTIONS, USE RESULT LATER
	$choose_starting_chr = mt_rand(0,1);
	
	//CREATE A NUMBER OF RECOMBINATION POSITIONS
	$rec_freq_current_chr = gauss_rec_freq($rec_freq, $rec_freq_sd); //Recombination frequency = $rec_freq = mean recombination events per base
	if ($rec_freq_current_chr < 0) {
		$rec_freq_current_chr = 0;
	}
	$rec_abs_count = round($genomeA_length * $rec_freq_current_chr);
	for($rec_points_count = 1; $rec_points_count <= $rec_abs_count;) {
		$rec_points[$rec_points_count - 1] = mt_rand(2, $genomeA_length - 1);
		$rec_points_count++;
	}
	$rec_points[] = 0;
	$rec_points[] = $genomeA_length;
	sort($rec_points);
	
	//DUPLICATE ARRAY WITH REC POSITIONS AND ADD STARTING CHROMOSOME INFO 
	$rec_point_and_starting_chr = $rec_points;
	array_unshift($rec_point_and_starting_chr, $choose_starting_chr);
/*
	echo "<pre>";
	print_r($rec_points);
	echo "</pre>";

	echo "<pre>";
	print_r($rec_point_and_starting_chr);
	echo "</pre>";
*/

	//STORE ALL GENOMES WITH THEIR RECOMBINATION POSITIONS AND STARTING CHR IN ARRAY TO DRAW IMAGE
	$array_for_image[$genome_count] = $rec_point_and_starting_chr;
	
	$ends_plus_rec_points_count = 2 + $rec_abs_count;
	//echo $rec_abs_count ."<br>";
	//echo "Recombinant genome number: ". $genome_count ."<br>";

	//CREATE RECOMBINANT GENOME USING RECOMBINATION POSITIONS
	$rec_genome = "";
	for($fragment_count = 1; $fragment_count < $ends_plus_rec_points_count;) {
		
		//RANDOMLY START WITH A GENOME
		if ($choose_starting_chr == 0) {
			if ($fragment_count & 1) { //simply a way of checking if number is odd or even
				$use_genome = $genomeA;
				$use_genome2 = "genomeA";
			} else {
				$use_genome = $genomeB;
				$use_genome2 = "genomeB";
			}
		} else {
			if ($fragment_count & 1) {
				$use_genome = $genomeB;
				$use_genome2 = "genomeB";
			} else {
				$use_genome = $genomeA;
				$use_genome2 = "genomeA";
			}
		}
		$fragment_length = $rec_points[$fragment_count] - $rec_points[$fragment_count - 1];
		$fragment = substr($use_genome, $rec_points[$fragment_count - 1], $fragment_length);
		//echo "info point ". $fragment_count ." a ". $rec_points[$fragment_count - 1] ." b ". $rec_points[$fragment_count] ." l ". $fragment_length ."<br>";
		
		//DETECT IF RECOMBINANT GENOME CARRIES MUTATION 1
		if ($use_genome2 == $mutation1_genome AND $rec_points[$fragment_count - 1] < $mutation1_position AND $mutation1_position <= $rec_points[$fragment_count]) {
			$carries_mutation1 = "YES";
		}
		
		//DETECT IF RECOMBINANT GENOME CARRIES MUTATION 2
		if ($use_genome2 == $mutation2_genome AND $rec_points[$fragment_count - 1] < $mutation2_position AND $mutation2_position <= $rec_points[$fragment_count]) {
			$carries_mutation2 = "YES";
		}
		
		//echo $fragment_count ."\t". $use_genome2 ."\t" . $rec_points[$fragment_count - 1] ."\t". $rec_points[$fragment_count] ."\t". $fragment_length ."<br>";
		$rec_genome .= $fragment;
		$fragment_count++;
	}
	
	//4 RECOMBINANT CHROMOSOME SELECTION MODES:
	if ($mutation_inherintance_mode == "r") {
		if (isset($carries_mutation1)) {
			//echo "mu ";
			//CREATE OUTPUT FILE WITH CURRENT MUTATION-CONTAINING RECOMBINANT GENOME
			$mut_rec_genome_fa = fopen("rec_genomes\genome_rec_". $genome_count ."_mut.fa","w");
			fwrite($mut_rec_genome_fa, $rec_genome);
			fclose($mut_rec_genome_fa);
			
			$genome_count++;
		}
	} else if ($mutation_inherintance_mode == "d") {
		if (isset($carries_mutation1)) {
			//echo "mu ";
			//CREATE OUTPUT FILE WITH CURRENT MUTATION-CONTAINING RECOMBINANT GENOME
			$mut_rec_genome_fa = fopen("rec_genomes\genome_rec_". $genome_count ."_mut.fa","w");
			fwrite($mut_rec_genome_fa, $rec_genome);
			fclose($mut_rec_genome_fa);
			
			$genome_count++;
		} else {
			//IN DOMINANT MUT-PHE MODE: WITH P = 0.5, CREATE OUTPUT FILE WITH CURRENT NON-MUTATION-CONTAINING RECOMBINANT GENOME
			$select_not_carries_mutation = mt_rand(1, 100);
			if ($select_not_carries_mutation > 50) { //P = 50 because in "dominant F2 population", picking mutant-phenotype inds will lead to picking 1/3 of non-mutation-containing chromosomes. Every 100 wt, take 50 mut, sum is 150, so 50/150 = 1/3
				$rec_genome_fa = fopen("rec_genomes\genome_rec_". $genome_count ."_notmut.fa","w");
				fwrite($rec_genome_fa, $rec_genome);
				fclose($rec_genome_fa);
				
				$genome_count++;
			}
			//echo "wt ";
		}		
	} else if ($mutation_inherintance_mode == "di") {
		if (!isset($carries_mutation1)) {
			//echo "wt ";
			//CREATE OUTPUT FILE WITH CURRENT MUTATION-CONTAINING RECOMBINANT GENOME
			$mut_rec_genome_fa = fopen("rec_genomes\genome_rec_". $genome_count ."_wt.fa","w");
			fwrite($mut_rec_genome_fa, $rec_genome);
			fclose($mut_rec_genome_fa);
			
			$genome_count++;
		}		
	} else if ($mutation_inherintance_mode == "dr") {
		if (isset($carries_mutation1) AND isset($carries_mutation2)) {
			//echo "mu ";
			//CREATE OUTPUT FILE WITH CURRENT MUTATION-CONTAINING RECOMBINANT GENOME
			$mut_rec_genome_fa = fopen("rec_genomes\genome_rec_". $genome_count ."_mut.fa","w");
			fwrite($mut_rec_genome_fa, $rec_genome);
			fclose($mut_rec_genome_fa);
			
			$genome_count++;
		}	
	}
	
	unset($rec_points);
	$carries_mutation1 = null; //bug: Some mut-containing chrs are selected!!! when I use unset($carries_mutation1) it happens the same.
	$carries_mutation2 = null;
	unset($select_not_carries_mutation);
	//Do this: $var = null to release memory
}

/*
echo "<pre>";
print_r($array_for_image);
echo "</pre>";
*/

//require_once('create_image.php');

?>

