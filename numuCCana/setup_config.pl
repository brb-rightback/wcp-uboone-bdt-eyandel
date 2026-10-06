#!/usr/bin/perl

my $new_config_string = $ARGV[0];
$new_config_string =~ s/configurations\///g;
print "New configuration is $new_config_string\n";
my @new_config = split(/\//, $new_config_string);

my $save = $ARGV[1] // 1;
if(${save} != 0){
  system("perl save_config.pl");
}

system("rm configurations/*.txt");

my $config_dir = "configurations";
foreach my $config (@new_config) {
  $config_dir = "${config_dir}/${config}";
  print("cp ${config_dir}/*.txt configurations/\n");
  system("cp ${config_dir}/*.txt configurations/");
}

$new_config_string =~ tr{/}{_};
system("rm config.txt");
open(FH, '>', "config.txt") or die $!;
print FH $new_config_string;
close(FH);
