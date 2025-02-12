use strict;
use warnings;
use Benchmark qw(:all);
use POSIX qw(strftime);
use FileHandle;
use Sys::Hostname;
use File::Basename qw( dirname );
use IO::Uncompress::UnXz qw(unxz $UnXzError) ;

use application 'polytope';
my $codename="collect_graphs_";

die "usage: $codename MPTOPCOM_INPUT MPTOPCOM_OUTPUT<n>\n" unless scalar(@ARGV)==2;

# Set some file names, read original points to later make sure that the order
# returned from polymake agrees with them.
my $input_file = $ARGV[0];
my $input = `cat $input_file`;
$input =~ s/\s//g;
my ($pts, $gp) = $input =~ m/(.*\]\])(\[\[.*)/;
$pts = new Matrix(eval($pts));
my $tdata=$ARGV[1];
my $output_dir = dirname(`realpath $tdata`);
print "Will write output to $output_dir\n";
my $root = `realpath $tdata`;
chomp $root;
$root =~ s/\.dat\.xz//;
my $output_file="$root.graphs";
my $output_classes = "$output_file.classes";
my $log_file="$output_file.log";
# Someone else is working on this file.
if(-e $output_file || -e $log_file || -e "$output_file.xz"){
   exit;
}

# mark start
my $host = hostname();
my $now = strftime "%a %b %e %H:%M:%S %Y", localtime();
my $LOG = FileHandle->new("> $log_file");
die "cannot write $log_file\n" unless defined($LOG);
$LOG->autoflush(1);
print $LOG "started $codename @ $host on $now\n";
print $LOG "reading $tdata\n";
# open my $TDATA,  "xzcat $tdata|" or die "cannot read $tdata";
my $TDATA = IO::Uncompress::UnXz->new($tdata)
    or die "IO::Uncompress::UnXz failed: $UnXzError\n";

print $LOG "writing $output_file\n";
my $OUT = FileHandle->new("> $output_file");
$OUT->autoflush(1);
die "cannot write $output_file\n" unless defined($OUT);


my $P = scale(simplex(3),2);
my $CP = cayley_polytope($P,$P);
my $cayley_points = $CP->LATTICE_POINTS;

# Check that the points are the same as given to mptopcom.
if($cayley_points->minor(All,~[5]) != $pts){
   die "Points do not agree with input points from mptopcom!";
}


my $minkowski = minkowski_sum($P, $P);
my $n_points = ($minkowski->LATTICE_POINTS)->rows();
my $points = $minkowski->LATTICE_POINTS;
my $ptsindices = new Map<Vector<Integer>, Int>();
for(my $i=0; $i<$points->rows(); $i++){
   $ptsindices->{$points->row($i)} = $i;
}
my @collected_graphs = ();

my %histogram = ();
my $classes = new Map<Int, GraphAdjacency>();
my $nclasses = 0;

#my $hash = new Array<Int>($size);

# start the clock
my $t0=Benchmark->new();

my $rayset1 = new Set<Vector<Int>>;


sub tropical_curve_from_subdivision{
   my($pts, $mc, $index_map) = @_;
   my $toblerone_cells = new Set<Set<Integer>>;

   foreach my $cell (@{$mc}){
      my $set1 = new Set<Vector<Rational>>;
      my $set2 = new Set<Vector<Rational>>;
      foreach my $vertex (@{$pts->minor($cell, All)}){
         if ($vertex->[4] == 1){$set1 += $vertex} else {$set2 += $vertex};}   
      if($set1->size()*$set2->size() == 6){
         my $sum_vert = new Set<Vector<Integer>>();
         foreach my $v1 (@$set1){
            foreach my $v2 (@$set2){
               my $s = $v1+$v2;
               $s->[0] = 1;
               $sum_vert += $s->slice(range(0,3));
            }
         }
         my $cell_indices = new Set<Integer>;
         foreach my $vertex (@{$sum_vert}){
            $cell_indices += $index_map->{$vertex};
         }
         $toblerone_cells += $cell_indices
      }
   }

   my $ordered_toblerones_cells = new Array<Set<Int>>($toblerone_cells);
   my $edges = new Set<Array<Int>>;
   for (my $i=0; $i<15; ++$i) {
      for (my $j=$i+1; $j<16; ++$j) {
         my $intersection = $ordered_toblerones_cells->[$i]*$ordered_toblerones_cells->[$j];
         if ($intersection->size() == 4) {$edges += [$i,$j];}
      }
   }
   return $edges;
}

my $id=0;
while(my $line=$TDATA->getline()) {
   chomp $line;
   $line =~ s/^.*(\{\{.*\}\}).*/$1/;
   $line =~ s/\{/\[/g; $line =~ s/\}/\]/g;
   my $T = new Array<Set>(eval($line));
   my $edges = tropical_curve_from_subdivision($cayley_points, $T, $ptsindices);
   my $g = graph_from_edges($edges);
   my $cycle_length = 0;
   foreach my $comp (@{$g->BICONNECTED_COMPONENTS}) {if ($comp->size>1) {$cycle_length = $comp->size; last;}}
   ++$histogram{$cycle_length};
#$hash->[$id] = $cycle_length; #which is faster?
   my $ch = canonical_hash($g->ADJACENCY);
   print $OUT "ID :$id; edges: $edges; canonical_hash: $ch\n";
   if($id%10000 == 0){
      my $ttmp = Benchmark->new();
      my $diff = timestr(timediff($ttmp,$t0));
      print $LOG "$id($nclasses; $diff)\n";
   }
   if(defined $classes->{$ch}){
      my $test = graph::isomorphic($classes->{$ch}, $g->ADJACENCY) ? 1:0;
      if($test == 0){
         die "Collision!\n";
# my @input = @{$classes->{$ch}};
# push @input, $g->ADJACENCY;
# my $newar = new Array<GraphAdjacency>(\@input);
# $classes->{$ch} = $newar;
# $nclasses++;
      }
   } else {
      $classes->{$ch} = $g->ADJACENCY;
      $nclasses++;
   }
   ++$id;
}

$TDATA->close();
my $t1=Benchmark->new();

close $OUT or die "cannot close $output_file\n";

my $td1=timediff($t1,$t0);

print $LOG "\n$id processed lines\n", "elapsed time for computation: ", timestr($td1), "\n";
print $LOG "histogram:\n";
foreach my $os (sort {$a <=> $b} keys(%histogram)) {
    print $LOG $os, ":", $histogram{$os}, " ";
}
save($classes, $output_classes);
# if($nclasses != $classes->size()){ die "Something went wrong when collecting classes!\n"}
print $LOG "\nCollected $nclasses graphs\n";
close $LOG or die "cannot close $log_file\n";
