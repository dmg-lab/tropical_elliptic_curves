use strict;
use warnings;
use Math::Trig;
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
   my($fn, $edges) = @_;
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
   foreach my $i (keys %xcoords){
      my $x = $xcoords{$i};
      my $y = $ycoords{$i};
      print $TIKZ "\\fill ($x, $y) circle (3pt);\n";
   }
   
   print $TIKZ "\\end{tikzpicture}\n";
   close $TIKZ or die "Cannot close tikz file $fn\n";
}

sub simplex_to_string {
   my($in) = @_;
   my %encoding = qw(0 A
   1 B
   2 C
   3 D
   4 E
   5 F
   6 G
   7 H
   8 I
   9 J
   10 K
   11 L
   12 M
   13 N
   14 O
   15 P
   16 Q
   17 R
   18 S
   19 T);
   my $result = join("", map($encoding{$_}, @$in));
   return $result;
}


my $id2edges = load("regular_unimodular.representative_graphs");
my $cl2ids = load("regular_unimodular.cl2ids");
my $classes_histogram = load("regular_unimodular.classes_histogram");
my $id2triang = load("regular_unimodular.representative_triangulations");

my $sum = 0;
my $counter = 0;
my $OUT = FileHandle->new("> inner.tex");
my ($mink,$minid, $maxk,$maxid) = (100000000,0,0,0);
if(!-d "graphs"){ `mkdir graphs`; }
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
      my $edges = $id2edges->{$id};
      my $triang = $id2triang->{$id};
      $triang = join(", ", map(simplex_to_string($_), @$triang));
      graph_to_tikz("graphs/$tikzfilename", $edges);
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
      print $OUT "Sample triangulation: \\tiny $triang \\normalsize \\\\\n";
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
