# QUAlloc

_QUAlloc_ is a globally applicable water use and allocation model that distributes the water demands to the water availability from different sources while considering the water quality requirements for the main water use sectors: domestic, irrigation, livestock, manufacturing and energy. _QUAlloc_ is a high-resolution global-scale model with 10 km (5 arcmin) spatial resolution and a flexible temporal resolution (i.e., daily, monthly or yearly).

_QUAlloc_ considers the changing of long-term water availability and water quality to meet the demand over time and evaluates how short-term variations in supply and demand may require adjustments to minimize any shortfalls per sector. Optionally, _QUAlloc_ considers the spatial distribution of water at any moment in time, aspects of the water infrastructure in terms of the network of sources and sinks, limiting abstraction potentials (i.e., pumping capacities), and water management considerations by prioritizing the order in which demands are met.

_QUAlloc_ is built up on the global hydrological model PCR-GLOBWB2 (Sutanudjaja et al., 2018; https://github.com/UU-Hydro/PCR-GLOBWB_model) and the global surface water quality model DynQual (Jones et al., 2023; https://github.com/SustainableWaterSystems/DYNQUAL).

Currently, _QUAlloc_ has a stand-alone modelling configuration with user-defined water demands, requiring hydrological and water quality input data. The model input files are:
- **water availability** (_obtained from PCR-GLOBWB2_): discharge, total runoff, direct runoff, interflow, baseflow, channel storage, groundwater storage, groundwater recharge and returns flow.
- **water quality** (_obtained from DynQual_): surface water temperature, biochemical oxygen demand, total dissolved solids, faecal coliforms.

The main output files from _QUAlloc_ are the water withdrawals and allocation volumes per water source (i.e., surface water, renewable groundwater and non-renewable groundwater) and per sector (i.e., domestic, irrigation, livestock, manufacturing and thermoelectric).

An self-contained example for the Rhine-Meuse basin can be found in https://doi.org/10.5281/zenodo.14511236.

**Contact (QUAlloc)**: Gabriel A. Cárdenas B. (g.a.cardenasbelleza@uu.nl)

**Contact (DynQual)**: Edward R. Jones (e.r.jones@uu.nl)

**Contact (PCR-GLOBWB2)**: Edwin H. Sutanudjaja (E.H.Sutanudjaja@uu.nl)

### References:
**DynQual**: Jones, E. R., Bierkens, M. F. P., Wanders, N., Sutanudjaja, E. H., van Beek, L. P. H., and van Vliet, M. T. H.: DynQual v1.0: a high-resolution global surface water quality model, Geosci. Model Dev., 16, 4481–4500, https://doi.org/10.5194/gmd-16-4481-2023, 2023.

**PCR-GLOBWB2**: Sutanudjaja, E. H., van Beek, R., Wanders, N., Wada, Y., Bosmans, J. H. C., Drost, N., van der Ent, R. J., de Graaf, I. E. M., Hoch, J. M., de Jong, K., Karssenberg, D., López López, P., Peßenteiner, S., Schmitz, O., Straatsma, M. W., Vannametee, E., Wisser, D., and Bierkens, M. F. P. (2018). PCR-GLOBWB 2: a 5 arcmin global hydrological and water resources model, Geosci. Model Dev., 11, 2429-2453, https://doi.org/10.5194/gmd-11-2429-2018.
