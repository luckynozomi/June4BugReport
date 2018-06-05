# June4iSegBugReport

## Bug Description
When using iSeg on a large `bc` value, the iSeg program will encounter an undefined behavior, which results in either
* The program outputs incorrectly, or
* The program outputs (almost) nothing.

## What iSeg Results is Affected
Technically, if after running function `iSeg::BenjaminiHochBerg()`, if the `SIG_CUTOFF` value (a.k.a. Final p-value cutoff) satisfies `SIG_CUTOFF <= 1e-10`, then the result will be affected.
Here is a table of `bc` value vs. typical `SIG_CUTOFF` value (the table is not precise, as it mixes different data)

| `bc` value | `SIG_CUTOFF`|
|:----------:|:-----------:|
|1.0|0.000275722|
|2.0|8.34952e-05|
|3.0|8.34952e-05|
|4.0|4.436e-06|
|5.0|1.55836e-07|
|6.0|5.49367e-11|
|7.0|5.48948e-13|
|8.0|1.54172e-17|
|9.0|4.32428e-20|
|9.6|4.23354e-22|
|10.0|5.82012e-24|

So iSeg runs with smaller `bc` value are not affected by this bug, while those with higher `bc` value will be affected.

## Technical Details
The bug lies in `iSeg::compSig` function called from `iSeg::output` function. Here is a snippet of the lines of `iSeg::compSig`

```
double iSeg::compSig(int start, int end, double* stat) {
	double pval_cut[12] = {1, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8};
	double pCutOff[10] = {0, 1.28, 2.326, 3.09, 3.719, 4.265, 4.753, 5.199, 5.612}; // test statistics cutoff for z-test
	double s1, s2, stat_t, stat_z, pval = 1, mean;
	int winlen, k, tmpWL;
	s1 = CS[end+1] - CS[start];
	s2 = CSS[end+1] - CSS[start];
	winlen = end - start + 1;
	if(CAL_TRUE_P == 0){
		NUM_TESTS++;
		stat_z = abs(s1) / sqrt(winlen);
		*stat = stat_z;
		return 1/stat_z;
	}

	mean = s1/(end-start+1);

	// use the cutoff to filter out non-significant regions
	if(abs(mean) <= BC * SD) {
		*stat = 0;
		return pval;
	}

	// calculate p-value using z-test

	stat_z = abs(s1) / (sqrt(winlen) * SD);

	if(COMP_Z == 1 && abs(stat_z) > pCutOff[int(-1*log10(SIG_CUTOFF))]){
		normal ndist(0, 1);
		pval = cdf(complement(ndist, abs(stat_z))) * 2; // two-tailed test
		NUM_TESTS++;
		*stat = stat_z;
		return pval;
	}
```

The expression `pCutOff[int(-1*log10(SIG_CUTOFF))]` failed to check whether `index = int(-1*log10(SIG_CUTOFF))` is in range [0, 10), as array `pCutOff` only has 10 elements.
Because C++ does not check whether or not an index is within range, when `index >= 10` the code will function in a different way than a human expects:
* When `10 <= index < 22`, because array `pval_cut` is defined earlier than `pCutOff`, it has a higher address in the stack, and thus `pCutOff[10] == pval_cut[0]` and so on so forth.
This makes condition `abs(stat_z) > pCutOff[int(-1*log10(SIG_CUTOFF))]` more easy to satisfy as intended, as those values in `pval_cut` are much smaller.
* When ` index >= 22`, the expression `pCutOff[index]` will access bits in memory that are not `double` formatted, and thus the resulting value will vary in different runs.

The discussion above and the relation between `bc` values and `SIG_CUTOFF` supports the observation that there is no error given in the entire pipeline when setting `bc` values to be 6.0, 7.0, 8.0, 8.5 and 9.0, but when bc = 9.9 or 10.0 (two values I tested), there is a random chance that the fusing right after running iSeg throws an error because the iSeg results are (almost) empty. 

## Suggested Fix

A quick fix is to extend the `pCutOff` array (but I don't know what values are those) to let it include more elements, and to warn the user when `index` is too big. The warning message will let the user know that the result is uncorrect because of it and he/she can adjust the array accordingly.

