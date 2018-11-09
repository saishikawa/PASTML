use strict;
use warnings;

my $taxa=1000;
open(NEWFILE, "> Constraint.txt") or die "$!";
print NEWFILE "(out,(";
for(my $i=1; $i <= $taxa; $i++){
  if($i < 1000){
    print NEWFILE "t$i,";
  } else {
    print NEWFILE "t$i));";
   }
}
close(NEWFILE);

