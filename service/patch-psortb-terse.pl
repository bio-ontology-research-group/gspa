# Preserve all tied best localizations. Upstream terse output selects the first
# Bayesian result, whose order can depend on Perl's randomized hash iteration.
use strict;
use warnings;
use Bio::Tools::PSort::Report::Formatter::terse;
my $path = $INC{'Bio/Tools/PSort/Report/Formatter/terse.pm'};
open(my $input, '<', $path) or die $!;
local $/; my $source = <$input>; close($input);
my $old = 'my ($res) = $rep->get_result("output", $_);';
my $new = 'my @results = $rep->get_result("output", $_); my $res = $results[0];';
$source =~ s/\Q$old\E/$new/ or die 'Unexpected PSORTb result formatter';
$old = 'my $loc = ($res->score >= $self->{cutoff}) ? $res->localization : "Unknown";';
$new = 'my $loc = ($res->score >= $self->{cutoff}) ? join("/", sort map { $_->localization } grep { abs($_->score - $res->score) < 1e-12 } @results) : "Unknown";';
$source =~ s/\Q$old\E/$new/ or die 'Unexpected PSORTb localization formatter';
open(my $output, '>', $path) or die $!; print $output $source; close($output);
