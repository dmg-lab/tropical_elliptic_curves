use strict;
use warnings;
use Math::Trig;
use Math::Complex;
use FileHandle;


use application "polytope";

sub format_number {
   my($n) = @_;
   if($n < 1000){
      return "$n";
   } else {
      my $tail = sprintf("%03d", $n%1000);
      my $head = format_number(($n-$tail)/1000);
      return "$head\\,$tail";
   }
}

sub graph_to_tikz {
   my($fn, $edges, $rn, $bn) = @_;
   my $g = graph_from_edges($edges);
   my $TIKZ = FileHandle->new("> $fn");
   print $TIKZ "\\begin{tikzpicture}[scale=\\tikzscale]\n";
   my $cycle_size = 0;
   my %xcoords = ();
   my %ycoords = ();
   for(my $i=0; $i<$g->N_NODES; $i++){
      $xcoords{$i} = 0;
      $ycoords{$i} = 0;
   }
   my @directed_cycle = ();
   foreach my $comp (@{$g->BICONNECTED_COMPONENTS}){
      if($comp->size > 2){
         my $start = $comp->front;
         push @directed_cycle, $start;
         my $oldnext = $start;
         my $next = ($g->ADJACENCY->adjacent_nodes($start) * $comp)->front;
         while($next != $start){
            push @directed_cycle, $next;
            my $candidates = ($g->ADJACENCY->adjacent_nodes($next) * $comp);
            $candidates -= $oldnext;
            $oldnext = $next;
            $next = $candidates->front;
         }
         # print join("->", @directed_cycle),"\n";
         
         $cycle_size = $comp->size();
         last;
      }
   }
   # Get cycle coords
   for(my $i=0; $i<$cycle_size; $i++){
      $xcoords{$directed_cycle[$i]} = sin(2*$i*pi/$cycle_size);
      $ycoords{$directed_cycle[$i]} = cos(2*$i*pi/$cycle_size);
   }
   my $i = 0;
   foreach my $node (@directed_cycle){
      my $done = new Set<Int>(@directed_cycle);
      my $neighbors = $g->ADJACENCY->adjacent_nodes($node) - $done;
      my @ordered_neighbors = @$neighbors;
      my $level = 2;
      my $levelstep = 0.9;
      while (@ordered_neighbors > 0){
         my @new_ordered_neighbors = ();
         my $step = 1.0/(scalar @ordered_neighbors + 1.0);
         my $k = 1.0;
         foreach my $on (@ordered_neighbors){
            $xcoords{$on} = $level*sin(2*($i-0.5+$k*$step)*pi/$cycle_size);
            $ycoords{$on} = $level*cos(2*($i-0.5+$k*$step)*pi/$cycle_size);
            $done += $on;
            my $onn = $g->ADJACENCY->adjacent_nodes($on) - $done;
            foreach my $nn (@$onn){
               push @new_ordered_neighbors, $nn;
            }
            $k++;
         }
         @ordered_neighbors = @new_ordered_neighbors;
         $level += $levelstep;
         $levelstep *= 0.9;
      }
      $i++;
   }
   my ($xmin,$xmax,$ymin,$ymax) = (0,0,0,0);
   foreach my $i (keys %xcoords){
      if($xcoords{$i} < $xmin){$xmin = $xcoords{$i};}
      if($ycoords{$i} < $ymin){$ymin = $ycoords{$i};}
      if($xcoords{$i} > $xmax){$xmax = $xcoords{$i};}
      if($ycoords{$i} > $ymax){$ymax = $ycoords{$i};}
   }
   if(($ymax-$ymin) > ($xmax-$xmin)){
      foreach my $i (keys %xcoords){
         my $tmp = $xcoords{$i};
         $xcoords{$i} = $ycoords{$i};
         $ycoords{$i} = $tmp;
      }
   }

   foreach my $e (@$edges){
      my $x0 = $xcoords{$e->[0]};
      my $y0 = $ycoords{$e->[0]};
      my $x1 = $xcoords{$e->[1]};
      my $y1 = $ycoords{$e->[1]};
      print $TIKZ "\\draw ($x0,$y0) -- ($x1,$y1);\n";
   }
   my $circle_size = (2*sqrt(16.0))/sqrt($cycle_size)+(1.0*(16-$cycle_size)/16.0);
   foreach my $i (keys %xcoords){
      my $x = $xcoords{$i};
      my $y = $ycoords{$i};
      my $color = $rn->contains($i) ? "red" : "blue";
      print $TIKZ "\\fill[$color] ($x, $y) circle (".$circle_size."pt);\n";
   }
   
   print $TIKZ "\\end{tikzpicture}\n";
   close $TIKZ or die "Cannot close tikz file $fn\n";
}

sub simplex_to_string {
   my($in) = @_;
   my %encoding = qw(0 A
   1 a
   2 B
   3 b
   4 C
   5 c
   6 D
   7 d
   8 E
   9 e
   10 F
   11 f
   12 G
   13 g
   14 H
   15 h
   16 I
   17 i
   18 J
   19 j);
   my $result = join("", map($encoding{$_}, @$in));
   return $result;
}

###############################################################################
# Function for computing graph from a triangulation.
sub tropical_curve_from_subdivision{
   my($pts, $mc, $index_map) = @_;
   my $toblerone_cells = new Set<Set<Integer>>;
   my $rednodes = new Set<Int>();
   my $bluenodes = new Set<Int>();
   my $redsimplices = new Set<Set<Int>>();
   my $bluesimplices = new Set<Set<Int>>();

   my $redtoblerones = new Set<Set<Int>>();

   foreach my $cell (@{$mc}){
      my $set1 = new Set<Vector<Rational>>;
      my $set2 = new Set<Vector<Rational>>;
      foreach my $vertex (@{$pts->minor($cell, All)}){
         if ($vertex->[4] == 1){$set1 += $vertex} else {$set2 += $vertex};}
      # There are only three possible cases for the convex hulls of $set1 and
      # $set2: a point, a line, or a triangle. In any case the elements of
      # $set1 and $set2 are the vertices of a maximal cell, hence these will be
      # vertices of the convex hulls as well. We only want those cases where
      # the Minkowski sum of these convex hulls is a toblerone, i.e. a prism
      # over a triangle. This happens for the case of a line and a triangle.
      # Thus the sizes of $set1 and $set2 must be 2 and 3, or the other way
      # around. Thus we can just test that the product of the sizes is 6, this
      # already excludes all other cases. To see that this is sufficient,
      # notice that $set1 and $set2 now form a disjoint subdivision of the
      # vertices of a 4-dim simplex.
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
         $toblerone_cells += $cell_indices;
         if($set1->size() == 3){
            $redtoblerones += $cell_indices;
            $redsimplices += $cell;
         } else {
            $bluesimplices += $cell;
         }
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
   for(my $i=0; $i<16; $i++){
      if($redtoblerones->contains($ordered_toblerones_cells->[$i])){
         $rednodes += $i;
      } else {
         $bluenodes += $i;
      }
   }
   return ($edges, $rednodes, $bluenodes, $redsimplices, $bluesimplices);
}

my $P = scale(simplex(3),2);
my $minkowski = minkowski_sum($P, $P);
my $CP = cayley_polytope($P,$P);
my $cayley_points = $CP->LATTICE_POINTS;
my $points = $minkowski->LATTICE_POINTS;
my $ptsindices = new Map<Vector<Integer>, Int>();
for(my $i=0; $i<$points->rows(); $i++){
   $ptsindices->{$points->row($i)} = $i;
}


my $id2edges = load("data/regular_unimodular.representative_graphs");
my $cl2ids = load("data/regular_unimodular.cl2ids");
my $classes_histogram = load("data/regular_unimodular.classes_histogram");
my $id2triang = load("data/regular_unimodular.representative_triangulations");

my $sum = 0;
my $counter = 0;
my $OUT = FileHandle->new("> latex/inner.tex");
my ($mink,$minid, $maxk,$maxid) = (100000000,0,0,0);
if(!-d "latex/graphs"){ `mkdir latex/graphs`; }
foreach my $cl (keys %$cl2ids){
   my $nids = $cl2ids->{$cl}->size();
   print "Running over cycle length $cl ($nids)\n";
   my $nisos = $cl2ids->{$cl}->size();
   $nisos = $nisos == 1 ? "$nisos class" : "$nisos classes";
   print $OUT "\\subsection{Cycle length $cl ($nisos)}\n\n";
   foreach my $id (@{$cl2ids->{$cl}}){
      my $class_members = $classes_histogram->{$id};
      if($class_members > $maxk){$maxid = $counter; $maxk=$class_members;}
      if($class_members < $mink){$minid = $counter; $mink=$class_members;}
      print "ID $id ($class_members)\n";
      $sum += $class_members;
      my $class_members_string = format_number($class_members);
      my $tikzfilename = $cl."_".$id."_graph.tikz";
      my $svgfilename = $cl."_".$id."_graph.svg";
      my $triang = $id2triang->{$id};
      my ($edges, $rn, $bn, $rs, $bs) = tropical_curve_from_subdivision($cayley_points, $triang, $ptsindices);
      $triang = join(", ", map{
            my $result = simplex_to_string($_);
         if($rs->contains($_)){
            $result = "\\textcolor{red}{$result}";
         } elsif($bs->contains($_)) {
            $result = "\\textcolor{blue}{$result}";
         }
         $result;
         }@$triang);
      graph_to_tikz("latex/graphs/$tikzfilename", $edges, $rn, $bn);
      my $idstring = format_number($id);
      print $OUT "\\fbox{\n";
      print $OUT "\\begin{minipage}{.45\\textwidth}\n";
      print $OUT "\\begin{tabular}{cc}\n";
      print $OUT "\\begin{minipage}{.5\\textwidth}\n";
      print $OUT "ID: $counter\\\\ CH: $idstring\\\\\n";
      print $OUT "\\small Members in class:\\\\ $class_members_string\\normalsize\n";
      print $OUT "\\end{minipage} & \n";
      print $OUT "\\begin{minipage}{.4\\textwidth}\n";
      print $OUT "\\begin{scaletikzpicturetowidth}{\\textwidth}\n";
      print $OUT "\\input{graphs/$tikzfilename}\n";
      print $OUT "\\end{scaletikzpicturetowidth}\n";
      print $OUT "\\end{minipage}\n";
      print $OUT "\\end{tabular}\\\\[.3cm]\n";
      print $OUT "{\\tiny\\begin{spacing}{1.1}\n$triang\n\\end{spacing} }\n";
      print $OUT "\\end{minipage}\n}\n\n";
      $counter++;
      # if($counter == 100){
      #    last;
      # }
   }
   # if($counter == 100){
   #    last;
   # }
}
# `find graphs/ -type f -print0 | xargs -0 sed -i 's/scale = 1/scale = \\\\tikzscale/g'`;
close $OUT or die "Cannot close inner.tex\n";
print "minid: $minid\n";
print "maxid: $maxid\n";
print "Total: $sum\n";
