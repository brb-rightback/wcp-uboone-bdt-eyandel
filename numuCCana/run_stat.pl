#!/usr/bin/perl

my $config = "NONE";
open my $file, '<', "config.txt" or die "Cannot open config.txt: $!";
my $config = <$file>;
chomp($config);
close $file;
print "Running MC stat for $config\n";

system("./bin/merge_hist -r0 -l0 -e2 > mcstat/0.log");

if($config =~ /val/ ){
  print "Running MC stat cor for $config\n";
  system("bin/stat_cov_matrix -r0 -h1&");
  system("bin/stat_pred_cov_matrix -r0 -h1");
  system("mv hist_rootfiles/run_pred_stat.root mcstat/");
  system("mv hist_rootfiles/run_data_stat.root mcstat/");
}

