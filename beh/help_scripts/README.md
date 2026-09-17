<h1>READ ME</h1>

<h2>Help scripts for data analysis</h2>

<h3>Get_staircase.m</h3>

+ 1<sup>st</sup> version by Florian Kasten: <code>tACSChallenge_ImportData.m</code>
+ 2<sup>nd</sup> version by Giuseppe Di Dona: <code>tACSChallenge_AnalyzeStaircase.m</code>
+ 3<sup>rd</sup> version by Wahriman Andrade de Araújo: <code>Get_staircase.m</code>

The new script adds the used Base Brightness (variable <code>a</code> in the tACS Challenge GUI) and Target Brightness (variable <code>b</code> in the tACS Challenge GUI) in the graph and table.

**Output example:**
<img width="1344" height="808" alt="sub-L19_S01_StaircaseResults_DATE_15_09_2026_12_45_31" src="https://github.com/user-attachments/assets/c4d4a512-7025-474f-b1a0-b8e41aef859a" />

<h3>swap_left_right_buttons_press.m</h3>

This script swaps the values of the variables <code>LeftButton</code> and <code>RightButton</code> inside the designated block files (e.g., <code>sub-L19_S01_B1_Sh_2026_03_03_14_32_06.tsv</code>). In addition, it creates a backup of the original files before the swap.

**Why?**

The button box used by lab 19 had a malfunctional right button. As a solution, during the data collection, we turned the button box in 180° allowing the participant to continue using their dominant hand (right hand) while pressing the assigned left button. In turn, during the data analysis, to avoid the miss assumption that the data was collected with only left-handed participants, we used a script to fix the participants handiness based on the demographic information provided by the them.

**Script's output**
```text
>> swap_left_right_buttons_press('Z:\scripts\tACSChallenge_Project\data\L19', 19, 1:9)

=== swap_left_right_buttons_press ===
data path : Z:\scripts\tACSChallenge_Project\data\L19
lab       : L19
subjects  : [1 2 3 4 5 6 7 8 9]
"values moved" = rows whose value was handed over to the other column
every row is swapped, whatever the value, so it must equal "rows"
------------------------------------------------------------------------------

sub-L19_S01
   B1  : Left -> 139253 values moved | Right -> 139253 values moved | 139253 rows | OK
   B2  : Left -> 450308 values moved | Right -> 450308 values moved | 450308 rows | OK
   B3  : Left -> 442443 values moved | Right -> 442443 values moved | 442443 rows | OK
   B4  : Left -> 455714 values moved | Right -> 455714 values moved | 455714 rows | OK
   B5  : Left -> 139253 values moved | Right -> 139253 values moved | 139253 rows | OK
   B6  : Left -> 446357 values moved | Right -> 446357 values moved | 446357 rows | OK
   B7  : Left -> 445947 values moved | Right -> 445947 values moved | 445947 rows | OK
   B8  : Left -> 449716 values moved | Right -> 449716 values moved | 449716 rows | OK
   B9  : Left -> 140132 values moved | Right -> 140132 values moved | 140132 rows | OK
   B10 : Left -> 440212 values moved | Right -> 440212 values moved | 440212 rows | OK
   B11 : Left -> 449470 values moved | Right -> 449470 values moved | 449470 rows | OK
   B12 : Left -> 450289 values moved | Right -> 450289 values moved | 450289 rows | OK
   B13 : Left -> 139586 values moved | Right -> 139586 values moved | 139586 rows | OK

sub-L19_S02
   B1  : Left -> 135982 values moved | Right -> 135982 values moved | 135982 rows | OK
   B2  : Left -> 454649 values moved | Right -> 454649 values moved | 454649 rows | OK
   B3  : Left -> 453338 values moved | Right -> 453338 values moved | 453338 rows | OK
   B4  : Left -> 440559 values moved | Right -> 440559 values moved | 440559 rows | OK
   B5  : Left -> 143082 values moved | Right -> 143082 values moved | 143082 rows | OK
   B6  : Left -> 444490 values moved | Right -> 444490 values moved | 444490 rows | OK
   B7  : Left -> 449734 values moved | Right -> 449734 values moved | 449734 rows | OK
   B8  : Left -> 433759 values moved | Right -> 433759 values moved | 433759 rows | OK
   B9  : Left -> 149718 values moved | Right -> 149718 values moved | 149718 rows | OK
   B10 : Left -> 458007 values moved | Right -> 458007 values moved | 458007 rows | OK
   B11 : Left -> 453239 values moved | Right -> 453239 values moved | 453239 rows | OK
   B12 : Left -> 439967 values moved | Right -> 439967 values moved | 439967 rows | OK
   B13 : Left -> 155938 values moved | Right -> 155938 values moved | 155938 rows | OK

sub-L19_S03
   B1  : Left -> 150859 values moved | Right -> 150859 values moved | 150859 rows | OK
   B2  : Left -> 443326 values moved | Right -> 443326 values moved | 443326 rows | OK
   B3  : Left -> 446439 values moved | Right -> 446439 values moved | 446439 rows | OK
   B4  : Left -> 457907 values moved | Right -> 457907 values moved | 457907 rows | OK
   B5  : Left -> 145452 values moved | Right -> 145452 values moved | 145452 rows | OK
   B6  : Left -> 442425 values moved | Right -> 442425 values moved | 442425 rows | OK
   B7  : Left -> 450289 values moved | Right -> 450289 values moved | 450289 rows | OK
   B8  : Left -> 442015 values moved | Right -> 442015 values moved | 442015 rows | OK
   B9  : Left -> 155036 values moved | Right -> 155036 values moved | 155036 rows | OK
   B10 : Left -> 456924 values moved | Right -> 456924 values moved | 456924 rows | OK
   B11 : Left -> 448077 values moved | Right -> 448077 values moved | 448077 rows | OK
   B12 : Left -> 444637 values moved | Right -> 444637 values moved | 444637 rows | OK
   B13 : Left -> 153807 values moved | Right -> 153807 values moved | 153807 rows | OK

sub-L19_S04
   B1  : Left -> 148489 values moved | Right -> 148489 values moved | 148489 rows | OK
   B2  : Left -> 445556 values moved | Right -> 445556 values moved | 445556 rows | OK
   B3  : Left -> 447596 values moved | Right -> 447596 values moved | 447596 rows | OK
   B4  : Left -> 458809 values moved | Right -> 458809 values moved | 458809 rows | OK
   B5  : Left -> 155364 values moved | Right -> 155364 values moved | 155364 rows | OK
   B6  : Left -> 452582 values moved | Right -> 452582 values moved | 452582 rows | OK
   B7  : Left -> 435379 values moved | Right -> 435379 values moved | 435379 rows | OK
   B8  : Left -> 441278 values moved | Right -> 441278 values moved | 441278 rows | OK
   B9  : Left -> 146762 values moved | Right -> 146762 values moved | 146762 rows | OK
   B10 : Left -> 441852 values moved | Right -> 441852 values moved | 441852 rows | OK
   B11 : Left -> 464134 values moved | Right -> 464134 values moved | 464134 rows | OK
   B12 : Left -> 457498 values moved | Right -> 457498 values moved | 457498 rows | OK
   B13 : Left -> 142584 values moved | Right -> 142584 values moved | 142584 rows | OK

sub-L19_S05
   B1  : Left -> 142755 values moved | Right -> 142755 values moved | 142755 rows | OK
   B2  : Left -> 455631 values moved | Right -> 455631 values moved | 455631 rows | OK
   B3  : Left -> 442361 values moved | Right -> 442361 values moved | 442361 rows | OK
   B4  : Left -> 425802 values moved | Right -> 425802 values moved | 425802 rows | OK
   B5  : Left -> 153808 values moved | Right -> 153808 values moved | 153808 rows | OK
   B6  : Left -> 444964 values moved | Right -> 444964 values moved | 444964 rows | OK
   B7  : Left -> 444719 values moved | Right -> 444719 values moved | 444719 rows | OK
   B8  : Left -> 449387 values moved | Right -> 449387 values moved | 449387 rows | OK
   B9  : Left -> 163720 values moved | Right -> 163720 values moved | 163720 rows | OK
   B10 : Left -> 442343 values moved | Right -> 442343 values moved | 442343 rows | OK
   B11 : Left -> 449388 values moved | Right -> 449388 values moved | 449388 rows | OK
   B12 : Left -> 454139 values moved | Right -> 454139 values moved | 454139 rows | OK
   B13 : Left -> 142339 values moved | Right -> 142339 values moved | 142339 rows | OK

sub-L19_S06
   B1  : Left -> 147172 values moved | Right -> 147172 values moved | 147172 rows | OK
   B2  : Left -> 444883 values moved | Right -> 444883 values moved | 444883 rows | OK
   B3  : Left -> 447176 values moved | Right -> 447176 values moved | 447176 rows | OK
   B4  : Left -> 451599 values moved | Right -> 451599 values moved | 451599 rows | OK
   B5  : Left -> 153152 values moved | Right -> 153152 values moved | 153152 rows | OK
   B6  : Left -> 452664 values moved | Right -> 452664 values moved | 452664 rows | OK
   B7  : Left -> 451436 values moved | Right -> 451436 values moved | 451436 rows | OK
   B8  : Left -> 454959 values moved | Right -> 454959 values moved | 454959 rows | OK
   B9  : Left -> 153562 values moved | Right -> 153562 values moved | 153562 rows | OK
   B10 : Left -> 450453 values moved | Right -> 450453 values moved | 450453 rows | OK
   B11 : Left -> 439639 values moved | Right -> 439639 values moved | 439639 rows | OK
   B12 : Left -> 459628 values moved | Right -> 459628 values moved | 459628 rows | OK
   B13 : Left -> 148647 values moved | Right -> 148647 values moved | 148647 rows | OK

sub-L19_S07
   B1  : Left -> 145620 values moved | Right -> 145620 values moved | 145620 rows | OK
   B2  : Left -> 451126 values moved | Right -> 451126 values moved | 451126 rows | OK
   B3  : Left -> 443262 values moved | Right -> 443262 values moved | 443262 rows | OK
   B4  : Left -> 445474 values moved | Right -> 445474 values moved | 445474 rows | OK
   B5  : Left -> 147343 values moved | Right -> 147343 values moved | 147343 rows | OK
   B6  : Left -> 449242 values moved | Right -> 449242 values moved | 449242 rows | OK
   B7  : Left -> 445628 values moved | Right -> 445628 values moved | 445628 rows | OK
   B8  : Left -> 446192 values moved | Right -> 446192 values moved | 446192 rows | OK
   B9  : Left -> 144387 values moved | Right -> 144387 values moved | 144387 rows | OK
   B10 : Left -> 440050 values moved | Right -> 440050 values moved | 440050 rows | OK
   B11 : Left -> 474373 values moved | Right -> 474373 values moved | 474373 rows | OK
   B12 : Left -> 444145 values moved | Right -> 444145 values moved | 444145 rows | OK
   B13 : Left -> 144059 values moved | Right -> 144059 values moved | 144059 rows | OK

sub-L19_S08
   B1  : Left -> 157090 values moved | Right -> 157090 values moved | 157090 rows | OK
   B2  : Left -> 433841 values moved | Right -> 433841 values moved | 433841 rows | OK
   B3  : Left -> 449715 values moved | Right -> 449715 values moved | 449715 rows | OK
   B4  : Left -> 438820 values moved | Right -> 438820 values moved | 438820 rows | OK
   B5  : Left -> 147580 values moved | Right -> 147580 values moved | 147580 rows | OK
   B6  : Left -> 448486 values moved | Right -> 448486 values moved | 448486 rows | OK
   B7  : Left -> 446193 values moved | Right -> 446193 values moved | 446193 rows | OK
   B8  : Left -> 441851 values moved | Right -> 441851 values moved | 441851 rows | OK
   B9  : Left -> 152906 values moved | Right -> 152906 values moved | 152906 rows | OK
   B10 : Left -> 446930 values moved | Right -> 446930 values moved | 446930 rows | OK
   B11 : Left -> 441360 values moved | Right -> 441360 values moved | 441360 rows | OK
   B12 : Left -> 451027 values moved | Right -> 451027 values moved | 451027 rows | OK
   B13 : Left -> 157330 values moved | Right -> 157330 values moved | 157330 rows | OK

sub-L19_S09
   B1  : Left -> 155943 values moved | Right -> 155943 values moved | 155943 rows | OK
   B2  : Left -> 464810 values moved | Right -> 464810 values moved | 464810 rows | OK
   B3  : Left -> 453830 values moved | Right -> 453830 values moved | 453830 rows | OK
   B4  : Left -> 452355 values moved | Right -> 452355 values moved | 452355 rows | OK
   B5  : Left -> 141362 values moved | Right -> 141362 values moved | 141362 rows | OK
   B6  : Left -> 450125 values moved | Right -> 450125 values moved | 450125 rows | OK
   B7  : Left -> 442424 values moved | Right -> 442424 values moved | 442424 rows | OK
   B8  : Left -> 456433 values moved | Right -> 456433 values moved | 456433 rows | OK
   B9  : Left -> 144960 values moved | Right -> 144960 values moved | 144960 rows | OK
   B10 : Left -> 448897 values moved | Right -> 448897 values moved | 448897 rows | OK
   B11 : Left -> 458399 values moved | Right -> 458399 values moved | 458399 rows | OK
   B12 : Left -> 442835 values moved | Right -> 442835 values moved | 442835 rows | OK
   B13 : Left -> 153316 values moved | Right -> 153316 values moved | 153316 rows | OK

------------------------------------------------------------------------------
files swapped: 117   files skipped: 0   mismatches: 0
all files    : Left -> 41642263 values moved | Right -> 41642263 values moved | 41642263 rows

==============================================================================
SANITY CHECK - last file modified: sub-L19_S09, block B13
Z:\scripts\tACSChallenge_Project\data\L19\sub-L19_S09\beh\sub-L19_S09_B13_Sh_2026_03_24_12_26_09.tsv
first rows with a button press (8095 press samples in total)
       row |   LeftButton  RightButton |   LeftButton  RightButton
           |       BEFORE       BEFORE |        AFTER        AFTER
--------------------------------------------------------------------
     13918 |            1            0 |            0            1
     13919 |            1            0 |            0            1
     13920 |            1            0 |            0            1
     13921 |            1            0 |            0            1
     13922 |            1            0 |            0            1
     13923 |            1            0 |            0            1
     13924 |            1            0 |            0            1
     13925 |            1            0 |            0            1
     13926 |            1            0 |            0            1
     13927 |            1            0 |            0            1
--------------------------------------------------------------------
whole file : Left -> 153316 values moved | Right -> 153316 values moved | 153316 rows
every row: LeftButton(after) == RightButton(before) and vice versa -> OK
==============================================================================
```
