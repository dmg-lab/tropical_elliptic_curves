use strict;
use warnings;
use Benchmark qw(:all);
use POSIX qw(strftime);
use FileHandle;
use Sys::Hostname;
use File::Basename qw( dirname );
use IO::Uncompress::UnXz qw(unxz $UnXzError) ;

use application 'polytope';
my $codename="extract_classes_";

die "usage: $codename MPTOPCOM_INPUT OUTPUT_PREFIX<n>\n" unless scalar(@ARGV)==2;

############################
# Helper methods
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
   if($toblerone_cells->size() != 16){
      die "Wrong number of toblerone cells!";
   }
   my $edges = new Set<Array<Int>>;
   for (my $i=0; $i<15; ++$i) {
      for (my $j=$i+1; $j<16; ++$j) {
         my $intersection = $ordered_toblerones_cells->[$i]*$ordered_toblerones_cells->[$j];
         if ($intersection->size() == 4) {$edges += [$i,$j];}
      }
   }
   return $edges;
}


sub old_tropical_curve_from_subdivision{
   my($pts, $mc, $index_map) = @_;
   my $subdivision = new fan::SubdivisionOfPoints(POINTS=>$pts, MAXIMAL_CELLS=>$mc);

   my $mixed_cells = new Set<Set<Integer>>;
   foreach my $cell (@{$subdivision->MAXIMAL_CELLS}){
      my $Q = new Polytope(POINTS=>$subdivision->POINTS->minor($cell, All));
      my $set1 = new Set<Vector<Rational>>;
      my $set2 = new Set<Vector<Rational>>;
      foreach my $vertex (@{$Q->VERTICES}){
         if ($vertex->[4] == 1){$set1 += $vertex} else {$set2 += $vertex};
      }
      my $Q1 = new Polytope(POINTS=>$set1);
      my $Q2 = new Polytope(POINTS=>$set2);
      my $minkowski_cell = project_full(minkowski_sum($Q1, $Q2));
      my $cell_indices = new Set<Integer>;
      foreach my $vertex (@{$minkowski_cell->VERTICES}){
         $cell_indices += $index_map->{$vertex};
      }
#$cells += $cell_indices;
      if ($cell_indices->size()==6){$mixed_cells += $cell_indices};
   }
   my $toblerones_cells = new Array<Set<Int>>($mixed_cells);
   if($toblerones_cells->size() != 16){
      die "Wrong number of toblerone cells!";
   }
   my $edges = new Set<Array<Int>>;
   for (my $i=0; $i<15; ++$i) {
      for (my $j=$i+1; $j<16; ++$j) {
         my $intersection = $toblerones_cells->[$i]*$toblerones_cells->[$j];
         if ($intersection->size() == 4) {
            $edges += [$i,$j];
         }
      }
   }
   return $edges;
}

sub get_cycle_length{
   my $edges = shift;
   my $g = graph_from_edges($edges);
   my $cycle_length = 0;
   foreach my $comp (@{$g->BICONNECTED_COMPONENTS}) {
      if ($comp->size>1) {
         $cycle_length = $comp->size; 
         last;
      }
   }
   return $cycle_length;
}


# Set some file names, read original points to later make sure that the order
# returned from polymake agrees with them.
my $input_file = $ARGV[0];
my $input = `cat $input_file`;
$input =~ s/\s//g;
my ($pts, $gp) = $input =~ m/(.*\]\])(\[\[.*)/;
$pts = new Matrix(eval($pts));
my $output_prefix=$ARGV[1];
my $output_dir = dirname(__FILE__);
print "Will write output to $output_dir\n";
my $log_file="$output_dir/$output_prefix.extract_classes.log";
my $output_file="$output_dir/$output_prefix.extract_classes";

if(-e $log_file){
   print "$log_file already exists! Exiting.\n";
   exit;
}


############################
# Start log
my $host = hostname();
my $now = strftime "%a %b %e %H:%M:%S %Y", localtime();
my $LOG = FileHandle->new("> $log_file");
die "cannot write $log_file\n" unless defined($LOG);
$LOG->autoflush(1);
print $LOG "started $codename @ $host on $now\n";
print $LOG "writing $output_file\n";
my $OUT = FileHandle->new("> $output_file");
$OUT->autoflush(1);
die "cannot write $output_file\n" unless defined($OUT);

############################
# Output containers
# Map canonical hash to a triangulation giving a graph with that canonical hash.
my $representative_triangulations = new Map<Int, Array<Set<Int>>>();
my $representative_triangulations_fn = "$output_dir/$output_prefix.representative_triangulations";
# For every canonical hash one graph with that canonical hash.
my $representative_graphs = new Map<Int, Set<Array<Int>>>();
my $representative_graphs_fn = "$output_dir/$output_prefix.representative_graphs";
# Map cycle length to the set of canonical hashes with this cycle length.
my $cycle_length_2_ids = new Map<Int, Set<Int>>();
my $cycle_length_2_ids_fn = "$output_dir/$output_prefix.cl2ids";
# Map canonical hash to number of times it appeared.
my $classes = new Map<Int, Int>();
my $classes_fn = "$output_dir/$output_prefix.classes_histogram";
print "Output in files:\n";
print $representative_triangulations_fn,"\n";
print $representative_graphs_fn,"\n";
print $cycle_length_2_ids_fn,"\n";


############################
# Geometric setup
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


# start the clock
my $t0=Benchmark->new();

############################
# Main loop
my @files = grep($_ =~ m/$output_prefix\.\d+\.graphs.xz/, `ls`);
map(chomp $_, @files);
my $total_counter = 0;
my $classes_counter = 0;
foreach my $file (@files){
   print $LOG "Reading file $file\n";
   my $TDATA = IO::Uncompress::UnXz->new($file)
      or die "IO::Uncompress::UnXz failed: $UnXzError\n";
   my $local_counter = 0;
   while(my $line = $TDATA->getline()){
      $total_counter++;
      # print $line;
      my($id, $edges, $ch) = $line =~ m/ID :(\d+); edges: (\{.*\}); canonical_hash: (\d+)/;
      $classes->{$ch} += 1;
      # print "ID: $id\n";
      $edges =~ s/\{/\[/;
      $edges =~ s/\}/\]/;
      $edges =~ s/<(\d+) (\d+)>/[$1,$2]/g;
      $edges =~ s/\] \[/\],\[/g;
      my $edges_str = $edges;
      $edges = new Set<Array<Int>>(eval $edges);
      # print "edges: $edges\n";
      # print "canonical_hash: $ch\n";
      if(!defined $representative_graphs->{$ch}){
         $classes_counter++;
         # print "New graph found $classes_counter!\n";
         $representative_graphs->{$ch} = $edges;
         my $triangfile = $file;
         $triangfile =~ s/\.graphs\.xz/\.dat\.xz/;
         my $tmpid = $id+1;
         my $triangline = `xzcat $triangfile | head -n $tmpid | tail -n1`;
         my ($triang) = $triangline =~ m/(\{\{.*\}\})/;
         $triang =~ s/\{/\[/g; $triang =~ s/\}/\]/g;
         my $triang_str = $triang;
         $triang = new Array<Set<Int>>(eval $triang);
         my $edges1 = tropical_curve_from_subdivision($cayley_points, $triang, $ptsindices);
         my $edges2 = old_tropical_curve_from_subdivision($cayley_points, $triang, $ptsindices);
         my $check = $edges1 == $edges2 && $edges == $edges1;
         if(!$check){ die "Methods disagree!" }
         $representative_triangulations->{$ch} = $triang;
         my $cl = get_cycle_length($edges);
         if(!defined $cycle_length_2_ids->{$cl}){
            $cycle_length_2_ids->{$cl} = new Set<Int>();
         }
         $cycle_length_2_ids->{$cl} += $ch;
         print $OUT "canonical_hash: $ch; triangulation: $triang_str; edges: $edges_str; cycle_length: $cl\n";
      }
      if($local_counter % 100000 == 0){
         print $LOG $local_counter, " ";
      }
      $local_counter++;
   }
   print $LOG "\n";
   $TDATA->close();
}

my $t1=Benchmark->new();
my $td1=timediff($t1,$t0);
print $LOG "Processed $total_counter lines\n";
print $LOG "Found $classes_counter classes of graphs.\n";
print $LOG "Elapsed time for computation: ", timestr($td1), "\n";
print $LOG "Saving output data!\n";
save($classes, $classes_fn);
save($cycle_length_2_ids, $cycle_length_2_ids_fn);
save($representative_triangulations, $representative_triangulations_fn);
save($representative_graphs, $representative_graphs_fn);
print $LOG "Done. Exiting.\n";

close $LOG or die "cannot close $log_file\n";
close $OUT or die "cannot close $output_file\n";

