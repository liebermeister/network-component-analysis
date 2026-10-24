function c = BEcolor()

c.TF      = [1 0 .7];
c.TFlight = 0.4*c.TF+0.6;

c.metabolite = [1 0 0];
c.metabolite_light = 0.4*c.metabolite+0.6;

c.flux = [.3 .3 .3];
c.flux_light = 0.4*c.flux+0.6;

c.protein = [1 0.7 0];
c.protein_light = 0.4*c.protein+0.6;

c.transcript = [0 0 1];
c.transcript_light = 0.4*c.transcript+0.6;