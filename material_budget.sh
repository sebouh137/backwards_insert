
phimin=-3.14159
phimax=3.14159
nphi=1


for xmlfile in beampipe_and_flanges_only.xml flange_-120cm_only.xml flange_-270cm_only.xml beampipe_flanges_and_ecal_endcap.xml; do
    for eta in -3.3 -3.35 -3.4 -3.45 -3.5 -3.55 -3.6 -3.65 -3.7 -3.75 -3.8 -3.85 -3.9 -3.95 -4.0 -4.05 -4.1 -4.15 -4.2 -4.25 -4.3 -4.35 -4.4 -4.45 -4.5 -4.55 -4.6 -4.65 -4.7 -4.75 -4.8 -4.85 -4.9; do
	python material_budget.py $xmlfile $eta $phimin $phimax $nphi &
    done
done
wait
echo \#done
