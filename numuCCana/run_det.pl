#!/usr/bin/perl

my $seed = 7;

my $config = "NONE";
open my $file, '<', "config.txt" or die "Cannot open config.txt: $!";
my $config = <$file>;
chomp($config);
close $file;
print "Running det sys for $config\n";

system("./bin/det_cov_matrix -r1 -h1 -s$seed&");
system("./bin/det_cov_matrix -r2 -h1 -s$seed&");
system("./bin/det_cov_matrix -r3 -h1 -s$seed&");
system("./bin/det_cov_matrix -r4 -h1 -s$seed&");
system("./bin/det_cov_matrix -r6 -h1 -s$seed&");
system("./bin/det_cov_matrix -r7 -h1 -s$seed&");
system("./bin/det_cov_matrix -r8 -h1 -s$seed&");
system("./bin/det_cov_matrix -r9 -h1 -s$seed&");
system("./bin/det_cov_matrix -r10 -h1 -s$seed&");
