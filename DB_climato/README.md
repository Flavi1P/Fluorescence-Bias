This folder contains all codes use to create the database of the NN (in project).<br> 
You can use the scripts in the order below.
<br><br>
**lov_hplc** : Shape the database of the HPLC measurements from the lov, create a csv<br><br>
**maredat_data** : Shape the database of the HPLC measurements from the maredat, create a csv<br><br>
**merge_datasets** : Combine the two database previously created, create a csv<br><br>
**par_climato** : Add climatology datas to the dataset + quality controle to make the definitiv dataset; create a csv <br>
What data did we chose ?<br>
<li>Only lov HPLC data</li>
<li>Profiles with more than 4 points</li>
<li>Profiles with the first measurement in the first 10m and the last one below the Zeu</li>
<li>Profiles with more than 500m depth</li>
<br><br>

