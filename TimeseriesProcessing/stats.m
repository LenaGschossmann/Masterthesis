% N-way ANOVA and Bonferroni post hoc
[f,t, st] = anovan(anova_val, {anova_grp1,anova_grp2},'model','interaction');
[c,m,h,g]=multcompare(st, 'Dimension',[2 1], 'CType', 'bonferroni');

% Test variables
x1=0;x2=0;x=0;

% 2-sample Kolmogorov-Smirnov test
[h,p,kst]=kstest2(x1,x2); % h=1: test rejecte H0 at 5%  significance level (H0 is that vectors are from same distribution)

% Pearson Correlation Coefficient
[r, p] = corrcoef(x(:,1), x(:,2));

% 1-way ANOVA
[f,t, st] = anova1(x(:,1), x(:,2));
[f,t, st] = anova1(anova_val, anova_grp);
[c,m,h,g]=multcompare(st, 'CType', 'bonferroni');
