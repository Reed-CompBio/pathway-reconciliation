```
python3 plot_ami.py corresponding-top-picks-7pathways-6overlap-withcollapsed.txt ../../out/
mv ../../out/combined_ami.pdf ../../out/combined_ami-6overlap.pdf

```
```
python3 plot_ami.py corresponding-top-picks-7pathways-4overlap-withcollapsed.txt ../../out/
mv ../../out/combined_ami.pdf ../../out/combined_ami-4overlap.pdf
```

```
python3 plot_ami_random.py netpath
python3 plot_ami_random_triangles2.py netpath

```

```
ls -d ../../networks/dbs/* | grep -v 'top-pick' | sed 's/..\/..\/networks\/dbs\///g' | awk '{print "python3 plot_ami_random.py "$1}' | bash
ls -d ../../networks/dbs/* | grep -v 'top-pick' | sed 's/..\/..\/networks\/dbs\///g' | awk '{print "python3 plot_ami_random_triangles2.py "$1}' | bash
```


## Withheld Results

```
python3 plot_dendrogram_leaveoneout.py ../../regression/predicted_withheld/dbs/*/Wnt_excluded* Wnt_excluded_ Wnt-excluded.png
```

To plot all:
```
 ls ../../regression/predicted_withheld/dbs/netpath/* | sed 's/excluded_.*/excluded_/g' | sed 's/.*netpath\///g' | sort -u | awk '{print "python3 plot_dendrogram_leaveoneout.py ../../regression/predicted_withheld/dbs/*/"$1"* "$1" figs-leaveoneout/"$1".png"}' | bash
```
