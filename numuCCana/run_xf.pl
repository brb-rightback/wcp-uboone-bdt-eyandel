#!/usr/bin/perl

my $config = "NONE";
open my $file, '<', "config.txt" or die "Cannot open config.txt: $!";
my $config = <$file>;
close $file;
chomp($config);
print "Running xf sys for $config\n";

print "Running flux sys";
system("./bin/xf_cov_matrix -r1 -h1&");
system("./bin/xf_cov_matrix -r2 -h1&");
system("./bin/xf_cov_matrix -r3 -h1&");
system("./bin/xf_cov_matrix -r4 -h1&");
system("./bin/xf_cov_matrix -r5 -h1&");
system("./bin/xf_cov_matrix -r6 -h1&");
system("./bin/xf_cov_matrix -r7 -h1&");
system("./bin/xf_cov_matrix -r8 -h1&");
system("./bin/xf_cov_matrix -r9 -h1&");
system("./bin/xf_cov_matrix -r10 -h1&");
system("./bin/xf_cov_matrix -r11 -h1&");
system("./bin/xf_cov_matrix -r12 -h1&");
system("./bin/xf_cov_matrix -r13 -h1&");

print "Running xs sys";
if($config =~ /unfold/ ){
  system("./bin/xs_cov_matrix -r17 -n0 -h1&");
}
 system("./bin/xf_cov_matrix -r17 -h1&");

print "Running reint sys";
system("./bin/xf_cov_matrix -r14 -h1&");
system("./bin/xf_cov_matrix -r15 -h1&");
system("./bin/xf_cov_matrix -r16 -h1&");
