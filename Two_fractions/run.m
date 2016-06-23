clear
swi=benthic_test.default_swi()
res=benthic_test.test_benthic(1,swi);
%benthic_test.plot_OMEN_BRNS(res, swi)
benthic_test.plot_column(res, false, swi)