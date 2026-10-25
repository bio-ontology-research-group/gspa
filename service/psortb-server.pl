use strict;
use warnings;
use HTTP::Daemon;
use HTTP::Response;
use File::Temp qw(tempdir);
use File::Path qw(make_path remove_tree);
use POSIX qw(setpgid);
$| = 1;
my $server = HTTP::Daemon->new(LocalAddr => '0.0.0.0', LocalPort => 8000, ReuseAddr => 1) or die $!;
print "PSORTb service listening\n";
while (my $client = $server->accept) {
    if (my $request = $client->get_request) {
        my $gram = $request->uri->query || '';
        if ($request->method ne 'POST' || $request->uri->path ne '/predict' || $gram !~ /^gram=(n|p|a)$/ || length($request->content) > 20971520) {
            $client->send_error(400);
        } else {
            my $mode = $1;
            my $dir = tempdir(CLEANUP => 1, DIR => '/tmp');
            open(my $input, '>', "$dir/input.faa") or die $!;
            print $input $request->content; close($input);
            my $pid = fork();
            if (!defined $pid) { die $!; }
            if ($pid == 0) {
                setpgid(0, 0);
                chdir($dir);
                open(STDOUT, '>', "$dir/output.tsv") or die $!;
                open(STDERR, '>', "$dir/error.log") or die $!;
                # PSORTb is serial. Four independent batches use the sidecar's
                # four CPU allocation without changing per-protein predictions.
                my @records = grep { /^>/ } split(/(?=^>)/m, $request->content);
                my $workers = @records < 4 ? scalar(@records) : 4;
                my @children;
                for my $part (0 .. $workers-1) {
                    make_path("$dir/part$part");
                    open(my $chunk, '>', "$dir/part$part/input.faa") or die $!;
                    for (my $i=$part; $i<@records; $i+=$workers) { print $chunk $records[$i]; }
                    close($chunk);
                    my $child = fork(); die $! unless defined $child;
                    if ($child == 0) {
                        chdir("$dir/part$part");
                        open(STDOUT, '>', 'output.tsv') or die $!;
                        open(STDERR, '>', 'error.log') or die $!;
                        exec('/usr/local/psortb/bin/psortb3.pl', '-'.$mode, '-o', 'terse', 'input.faa');
                        exit(1);
                    }
                    push @children, $child;
                }
                my $failed = 0;
                for my $child (@children) { waitpid($child,0); $failed=1 if $?; }
                for my $part (0 .. $workers-1) {
                    open(my $error, '<', "$dir/part$part/error.log"); print STDERR $_ while <$error>; close($error);
                    open(my $rows, '<', "$dir/part$part/output.tsv") or die $!;
                    my $header=<$rows>; print $header if $part==0 && defined $header;
                    print $_ while <$rows>; close($rows);
                }
                exit($failed);

            }
            my $timedout = 0;
            local $SIG{ALRM} = sub { $timedout = 1; kill 'KILL', -$pid; };
            alarm(3500); waitpid($pid, 0); my $status = $?; alarm(0);
            open(my $out, '<', "$dir/output.tsv") or die $!;
            local $/; my $body = <$out>; close($out);
            if ($status != 0 || $timedout) {
                open(my $err, '<', "$dir/error.log"); my $error = <$err>; close($err);
                print STDERR $error || 'PSORTb failed';
                $client->send_error(500, 'PSORTb did not complete');
            } else {
                $client->send_response(HTTP::Response->new(200, 'OK', ['Content-Type' => 'text/tab-separated-values'], $body));
            }
            remove_tree($dir);
        }
    }
    $client->close;
}
