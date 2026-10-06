#!/usr/bin/perl

my $old_config = "NONE";
open my $file, '<', "config.txt" or die "Cannot open config.txt: $!";
my $old_config = <$file>;
chomp($old_config);
print "The old config is $old_config\n";
close $file;

my $counter = 0;
while($counter<100){
  my $directory_path = "hist_rootfiles_$old_config\_$counter";
  if(-d $directory_path){
   $counter++; 
  }else{
    last;
  }
}
if ($counter==100) {
  print "Too many hist_rootfiles";
  exit 1;
}

print "Moving old files to hist_rootfiles_$old_config\_$counter\n";
system("mv hist_rootfiles hist_rootfiles_$old_config\_$counter");
system("mv mcstat hist_rootfiles_$old_config\_$counter/");
system("mv merge_xs.root hist_rootfiles_$old_config\_$counter/");
system("mv merge.root hist_rootfiles_$old_config\_$counter/");
system("mv file_collapsed_covariance_matrix.root hist_rootfiles_$old_config\_$counter/");
system("mkdir hist_rootfiles_$old_config\_$counter/plots/");
system("mv *.png hist_rootfiles_$old_config\_$counter/plots/");
system("mv *.pdf hist_rootfiles_$old_config\_$counter/plots/");
system("mv *.root hist_rootfiles_$old_config\_$counter/plots/");
system("mkdir hist_rootfiles_$old_config\_$counter/configuration/");
system("cp configurations/*.txt hist_rootfiles_$old_config\_$counter/configuration/");

system("mkdir hist_rootfiles");
system("mkdir hist_rootfiles/DetVar/");
system("mkdir hist_rootfiles/XsFlux/");
system("mkdir mcstat");

