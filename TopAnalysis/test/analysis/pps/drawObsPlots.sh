python test/analysis/pps/drawPvalCurve.py PPzX:0:ppvx_2017_unblind_signed_obs/optim_46 PPzmmX:0:ppvx_2017_unblind_signed_obs/optim_46 PPzeeX:0:ppvx_2017_unblind_signed_obs/optim_46
python test/analysis/pps/showFitShapes.py ppvx_2017_unblind_signed_obs/optim_46 0 z
python test/analysis/pps/mergeCrossingAngleDatacards.py ppvx_2017_unblind_signed_obs/optim_46
python test/analysis/pps/showFitShapes.py ppvx_2017_unblind_signed_obs/optim_46/inclusive 0 z
mv *.{png,pdf} /eos/user/p/psilva/www/ExclusiveAna_2017_unblind/obs/highpt_signed/
mv limits*root /eos/user/p/psilva/www/ExclusiveAna_2017_unblind/obs/highpt_signed/

python test/analysis/pps/drawPvalCurve.py PPzX:0:ppvx_2017_unblind_signed_obs/optim_10 PPzmmX:0:ppvx_2017_unblind_signed_obs/optim_10 PPzeeX:0:ppvx_2017_unblind_signed_obs/optim_10 PPgX:0:ppvx_2017_unblind_signed_obs/optim_2
python test/analysis/pps/showFitShapes.py ppvx_2017_unblind_signed_obs/optim_10 0 z
python test/analysis/pps/showFitShapes.py ppvx_2017_unblind_signed_obs/optim_2 0 g
python test/analysis/pps/mergeCrossingAngleDatacards.py ppvx_2017_unblind_signed_obs/optim_10
python test/analysis/pps/showFitShapes.py ppvx_2017_unblind_signed_obs/optim_10/inclusive 0 z
python test/analysis/pps/mergeCrossingAngleDatacards.py ppvx_2017_unblind_signed_obs/optim_2
python test/analysis/pps/showFitShapes.py ppvx_2017_unblind_signed_obs/optim_2/inclusive 0 g
mv *.{png,pdf} /eos/user/p/psilva/www/ExclusiveAna_2017_unblind/obs/signed/
mv limits*root /eos/user/p/psilva/www/ExclusiveAna_2017_unblind/obs/signed/


python test/analysis/pps/drawPvalCurve.py PPzX:0:ppvx_2017_unblind_obs/optim_46 PPzmmX:0:ppvx_2017_unblind_obs/optim_46 PPzeeX:0:ppvx_2017_unblind_obs/optim_46
python test/analysis/pps/showFitShapes.py ppvx_2017_unblind_obs/optim_46 0 z
python test/analysis/pps/mergeCrossingAngleDatacards.py ppvx_2017_unblind_obs/optim_46
python test/analysis/pps/showFitShapes.py ppvx_2017_unblind_obs/optim_46/inclusive 0 z
mv *.{png,pdf} /eos/user/p/psilva/www/ExclusiveAna_2017_unblind/obs/highpt/
mv limits*root /eos/user/p/psilva/www/ExclusiveAna_2017_unblind/obs/highpt/

python test/analysis/pps/drawPvalCurve.py PPzX:0:ppvx_2017_unblind_obs/optim_10 PPzmmX:0:ppvx_2017_unblind_obs/optim_10 PPzeeX:0:ppvx_2017_unblind_obs/optim_10 PPgX:0:ppvx_2017_unblind_obs/optim_2
python test/analysis/pps/showFitShapes.py ppvx_2017_unblind_obs/optim_10 0 z
python test/analysis/pps/showFitShapes.py ppvx_2017_unblind_obs/optim_2 0 g
python test/analysis/pps/mergeCrossingAngleDatacards.py ppvx_2017_unblind_obs/optim_10
python test/analysis/pps/showFitShapes.py ppvx_2017_unblind_obs/optim_10/inclusive 0 z
python test/analysis/pps/mergeCrossingAngleDatacards.py ppvx_2017_unblind_obs/optim_2
python test/analysis/pps/showFitShapes.py ppvx_2017_unblind_obs/optim_2/inclusive 0 g
mv *.{png,pdf} /eos/user/p/psilva/www/ExclusiveAna_2017_unblind/obs/
mv limits*root /eos/user/p/psilva/www/ExclusiveAna_2017_unblind/obs/


python test/analysis/pps/compareLimits.py \
    'baseline':/eos/user/p/psilva/www/ExclusiveAna_2017_unblind/obs/limits_PPzX_0.root \
    'high p_{T}':/eos/user/p/psilva/www/ExclusiveAna_2017_unblind/obs/highpt/limits_PPzX_0.root \
    'y(Z)-signed':/eos/user/p/psilva/www/ExclusiveAna_2017_unblind/obs/signed/limits_PPzX_0.root \
    'high p_{T},y(Z)-signed':/eos/user/p/psilva/www/ExclusiveAna_2017_unblind/obs/highpt_signed/limits_PPzX_0.root
mv limit*.{png,pdf} /eos/user/p/psilva/www/ExclusiveAna_2017_unblind/obs
