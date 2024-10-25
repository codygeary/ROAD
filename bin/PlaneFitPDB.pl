use strict;
use warnings FATAL => qw ( all );
use Data::Dumper;
use Math::Trig ':radial';
use Math::Trig;

# >>>>>>>>>>>>>>>>>> RUN PROGRAM <<<<<<<<<<<<<<<<<<<<

my ( $file, $line, @cols, $ATOM, @ATOM );

( $file ) = @ARGV;

if ( not open FILE, "< $file" ) {
    die "file not found!";
}

$ATOM = 0;
@ATOM = ( );
while ( $line = <FILE> ) {
    @cols = ( );
    @cols = split(/\s+/, $line);
    if ( defined $cols[0] ) {
        if ( $cols[0] eq "ATOM") {
            push @ATOM, {
                "a" => substr("$line", 13, 3),  #atom type
                "n" => substr("$line", 19, 1),  #residue type
                "i" => int(substr("$line", 22, 4)),  #residue number
                "h" => substr("$line", 21, 1), #chain name
                "x" => substr("$line", 30, 8),
                "y" => substr("$line", 38, 8),
                "z" => substr("$line", 46, 8),
                "j" => "N",    #set each atom to N, which means it might be an O3' position
            };
        }
        if ( $cols[0] eq "HETATM") {
            push @ATOM, {
                "a" => substr("$line", 13, 3),  #atom type
                "n" => substr("$line", 17, 1),  #residue type  ##for HETATM the residue type always gets printed 2 positions earlier in swisspdb!
                "i" => int(substr("$line", 22, 4)),  #residue number
                "h" => substr("$line", 21, 1), #chain name
                "x" => substr("$line", 30, 8),
                "y" => substr("$line", 38, 8),
                "z" => substr("$line", 46, 8),
                "j" => "N",    #set each atom to N, which means it might be an O3' position
            };
        }
    }
}

#print Dumper ( @ATOM );
#exit;
# PRINT PDB FILE

my $i = 0;
my $x = 1;
my $y = 1;
my $z = 1;
my $c = $ATOM[0]->{i};
my $W = 1;
my $AT = "";
my $t = 1;
my $residue=1;
my $thresh=0;
my $distance = 100;
my $cur_chain = "";
my $verbose = 1; #toggle verbose mode  1=verbose, 0=quiet

my $c_one_x=0; my $c_one_y=0; my $c_one_z=0;  #coords of all the C1-primes
my $c_four_x=0; my $c_four_y=0; my $c_four_z=0;  #coords of all the C4-primes
my $n_x=0; my $n_y=0; my $n_z=0;  #coords of all the N1 or N9

#------
#Calculate these from the PDB:
my $sum_zz = 0; my $n_sum = 0;  #number of atoms
my $com_x = 0;
my $com_y = 0;
my $com_z = 0;  #these are sums of the X,Y,Z coords that will set the translation vector for centering the PDB
#------

# first find P and C3* of first nucleotide
# store P and C3* xyz
# then look at next nucleotide

print "REMARK ---File Generated by PLaneFitPDB.pl - written by Cody Geary \n";

# Count residues

foreach $ATOM ( @ATOM ) {      #counts assuming the residues are all numbered differently
    if ( $ATOM->{i} != $c ) {
        $residue++;
        $c = $ATOM->{i};
    }
}
if ($verbose ==1) {print "REMARK "; print $residue; print " Residues in the PDB file \n";}

#$c = $ATOM[0]->{i};   #start at the first residue in the file
$c = 1;  #start at the residue numbered 1

$cur_chain = $ATOM[0]->{h};
if ($verbose ==1) {print "REMARK starting on chain: "; print $cur_chain; print "\n";}


foreach $ATOM ( @ATOM ) {  #First, fully calculate out the center-of-mass.
    $x = $ATOM->{x};
    $y = $ATOM->{y};
    $z = $ATOM->{z};

    $n_sum += 1;  #count up how many atoms total

    $com_x += $x;
    $com_y += $y;
    $com_z += $z;  #these are sums of the X,Y,Z coords that will set the translation vector for centering the PDB
}

$com_x = $com_x/$n_sum;  #now that they are counted up, divide them by the total number of atoms to get the center
$com_y = $com_y/$n_sum;
$com_z = $com_z/$n_sum;

print ("REMARK - center of mass : $com_x $com_y $com_z \n");
#center the coordinates and rotate_x to minimize, and then rotate_y to minimize sum(Z^2)/N

$sum_zz=0;
my $best_sum_zz = 0;  #set the initial sum very large
my $rotate_x_angle = 0;
my $rotate_y_angle = 0;

my $best_x = 0;
my $best_y = 0;
my $current_angle =0;


for (my $k = 0; $k<12; $k +=1){
    for (my $l = 0; $l<12; $l+=1){

        my $rand_x = ($k/12)*3.141592654;
        my $rand_y = ($l/12)*3.141592654;


        foreach $ATOM ( @ATOM ) {
            if ( $ATOM->{a} eq "P  ") {   #look only at phosphorous to save time
                #with the center-of-mass calculated, we can subtract it from each x,y,z for the sums
                $x = $ATOM->{x}-$com_x;
                $y = $ATOM->{y}-$com_y;
                $z = $ATOM->{z}-$com_z;

                my @trial = ();
                $trial[0]=$x;
                $trial[1]=$y;
                $trial[2]=$z;
                $trial[3]=1;

                @trial = &rotate_x ("@trial","$rand_x");
                @trial = &rotate_y ("@trial","$rand_y");
                $z = $trial[2];

                $sum_zz += abs($z);
            }
        }


        if ($best_sum_zz == 0 ){$best_sum_zz = $sum_zz;}

        if ($sum_zz < $best_sum_zz){
            $best_x = $rand_x;
            $best_y = $rand_y;
            $best_sum_zz = $sum_zz;
            print ("REMARK $k- X-angle = $rand_x  Y-angle = $rand_y  sum $sum_zz  best $best_sum_zz \n");
        }
        $sum_zz =0;
    }
}

my $very_best_x = $best_x;
my $very_best_y = $best_y;

for (my $k = 0; $k<12; $k +=1){
    for (my $l = 0; $l<12; $l+=1){

        my $rand_x = (($k-6)/6)*(1/12)*3.141592654+$best_x;
        my $rand_y = (($l-6)/6)*(1/12)*3.141592654+$best_y;


        foreach $ATOM ( @ATOM ) {
            if ( $ATOM->{a} eq "P  ") {   #look only at phosphorous to save time
                #with the center-of-mass calculated, we can subtract it from each x,y,z for the sums
                $x = $ATOM->{x}-$com_x;
                $y = $ATOM->{y}-$com_y;
                $z = $ATOM->{z}-$com_z;

                my @trial = ();
                $trial[0]=$x;
                $trial[1]=$y;
                $trial[2]=$z;
                $trial[3]=1;

                @trial = &rotate_x ("@trial","$rand_x");
                @trial = &rotate_y ("@trial","$rand_y");
                $z = $trial[2];

                $sum_zz += abs($z);
            }
        }


        if ($best_sum_zz == 0 ){$best_sum_zz = $sum_zz;}

        if ($sum_zz < $best_sum_zz){
            $very_best_x = $rand_x;
            $very_best_y = $rand_y;
            $best_sum_zz = $sum_zz;
            print ("REMARK $k- X-angle = $rand_x  Y-angle = $rand_y  sum $sum_zz  best $best_sum_zz \n");
        }
        $sum_zz =0;
    }
}

$best_x = $very_best_x;
$best_y = $very_best_y;  #get rid of my dumb temporary holders


for (my $k = 0; $k<12; $k +=1){
    for (my $l = 0; $l<12; $l+=1){

        my $rand_x = (($k-6)/6)*(1/144)*3.141592654+$best_x;
        my $rand_y = (($l-6)/6)*(1/144)*3.141592654+$best_y;


        foreach $ATOM ( @ATOM ) {
            if ( $ATOM->{a} eq "P  ") {   #look only at phosphorous to save time
                #with the center-of-mass calculated, we can subtract it from each x,y,z for the sums
                $x = $ATOM->{x}-$com_x;
                $y = $ATOM->{y}-$com_y;
                $z = $ATOM->{z}-$com_z;

                my @trial = ();
                $trial[0]=$x;
                $trial[1]=$y;
                $trial[2]=$z;
                $trial[3]=1;

                @trial = &rotate_x ("@trial","$rand_x");
                @trial = &rotate_y ("@trial","$rand_y");
                $z = $trial[2];

                $sum_zz += abs($z);
            }
        }

        if ($best_sum_zz == 0 ){$best_sum_zz = $sum_zz;}

        if ($sum_zz < $best_sum_zz){
            $very_best_x = $rand_x;
            $very_best_y = $rand_y;
            $best_sum_zz = $sum_zz;
            print ("REMARK $k- X-angle = $rand_x  Y-angle = $rand_y  sum $sum_zz  best $best_sum_zz \n");
        }
        $sum_zz =0;
    }
}

$best_x = $very_best_x;
$best_y = $very_best_y;  #get rid of my dumb temporary holders

print ("REMARK  best x=$best_x and best y=$best_y  and my best sumzz was $best_sum_zz \n");

foreach $ATOM ( @ATOM ) {
    if ( $ATOM->{j} eq "N" ) {
        printf "%-6s", "ATOM";
        printf "%5s", $t;
        print "  ";
        print $ATOM->{a};
        print "   ";
        print $ATOM->{n};
        print " ";
        ##print $ATOM->{h};  #Retain chain names
        print "A";           #or, force all chains to be A
        printf "%4s", $ATOM->{i};
        print "    ";

        my $coord_x = $ATOM->{x}-$com_x;
        my $coord_y = $ATOM->{y}-$com_y;
        my $coord_z = $ATOM->{z}-$com_z;

        ## Here we translate and rotate our final data set:

        my @final = ();   #Current Atom
        $final[0]=$coord_x;
        $final[1]=$coord_y;
        $final[2]=$coord_z;
        $final[3]=1;

        @final = &rotate_x ("@final",$best_x);  #first rotation
        @final = &rotate_y ("@final",$best_y); #second rotation

        my $roundedx =0; my $roundedy=0; my $roundedz=0;

        $roundedx = sprintf("%4.3f", $final[0]);
        $roundedy = sprintf("%4.3f", $final[1]);
        $roundedz = sprintf("%4.3f", $final[2]);

        printf("%8s", $roundedx);
        printf("%8s", $roundedy);
        printf("%8s", $roundedz);

                print "  1.00100.00   ";
        print "\n";
        $t++;

        $ATOM->{j} = "X";  #Mark the position as found
    }
}

print "END 1\n";


#-------------------------------subroutines

sub translate_matrix {
    my @r_tr1 = split(' ',$_[0]);
    my @r_tr2 = split(' ',$_[1]);
    my $tx=$r_tr2[0];
    my $ty=$r_tr2[1];
    my $tz=$r_tr2[2];
    my @Translate=(
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        $tx, $ty, $tz, 1);

    transform ("@r_tr1","@Translate");
}

sub rotate_x {  #point followed by theta in rad
    my @r_tr1 = split(' ',$_[0]);
    my $theta = $_[1];
    my @Translate=(
        1, 0, 0, 0,
        0, cos($theta), -sin($theta), 0,
        0, sin($theta), cos($theta), 0,
        0, 0, 0, 1);
    transform ("@r_tr1","@Translate");
}

sub rotate_y {  #point followed by theta in rad
    my @r_tr1 = split(' ',$_[0]);
    my $theta = $_[1];
    my @Translate=(
        cos($theta), 0, sin($theta), 0,
        0, 1, 0, 0,
        -sin($theta), 0, cos($theta), 0,
        0, 0, 0, 1);
    transform ("@r_tr1","@Translate");
}

sub rotate_z {  #point followed by theta in rad
    my @r_tr1 = split(' ',$_[0]);
    my $theta = $_[1];
    my @Translate=(
        cos($theta), -sin($theta), 0, 0,
        sin($theta), cos($theta), 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1);
    transform ("@r_tr1","@Translate");
}

sub transform {   #point (x,y,z,1) multiplied by the transform_matrix (16 elements long)
    my @m1 = split(' ',$_[0]);
    my @m2 = split(' ',$_[1]);
    my $result = [];

    my $xout=$m1[0]*$m2[0] + $m1[1]*$m2[4] + $m1[2]*$m2[8] + $m1[3]*$m2[12];
    my $yout=$m1[0]*$m2[1] + $m1[1]*$m2[5] + $m1[2]*$m2[9] + $m1[3]*$m2[13];
    my $zout=$m1[0]*$m2[2] + $m1[1]*$m2[6] + $m1[2]*$m2[10] + $m1[3]*$m2[14];
    my $vout=$m1[0]*$m2[3] + $m1[1]*$m2[7] + $m1[2]*$m2[11] + $m1[3]*$m2[15];  #we don't use the last digit, so no need to calculate it
    ($xout, $yout, $zout, $vout);
}
