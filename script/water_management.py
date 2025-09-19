# water management module of the QUAlloc model

# modules
import sys
import logging
import pcraster as pcr

from copy            import deepcopy
from basic_functions import pcr_return_val_div_zero, sum_list, pcr_get_statistics, max_dicts
from model_time      import match_date_by_julian_number, get_weights_from_dates
from allocation      import get_key, get_zonal_fraction, get_zonal_total, \
                            obtain_allocation_ratio, \
                            allocate_demand_to_availability_with_options, \
                            allocate_demand_to_withdrawals
from water_quality   import water_quality
#from estimate_waterdepth import estimate_waterdepth_from_discharge

# global attributes
# set the logger
logger = logging.getLogger(__name__)

# type of None, compatible with python 2.6
NoneType = type(None)

# small number to avoid zero divisions in PCRaster
very_small_number = 1.00e-12

# missing_id
mid = -1

# missing value identifier for water management options, including
# long-term water availability; missig values are then generated that are
# then treated in a default value to allocate water demand to resources and
# update long-term water availability etc.
water_management_missing_value = -9.99

# path out for debugging
path = '/scratch/carde003/qualloc/_debug'
verbose = False

########
# TODO #
########
critical_improvements= str.join('\n\t',\
             ( \
              '',\
              ))

development= str.join('\n\t',\
             ( \
              '', \
              ))

print ('\nDevelopmens for main module:')

if len(critical_improvements) > 0:
    print('Critical improvements: \n%s' % \
          critical_improvements)

if len(development) > 0:
    print ('Ongoing: \n%s' % development)

if len(critical_improvements) > 0:
    sys.exit()


#############
# functions #
#############

def estimate_waterdepth_from_discharge(discharge, \
                                       waterdepth, \
                                       mannings_n, \
                                       channel_width, \
                                       channel_gradient, \
                                       beta = 0.60, \
                                       max_iterations = 25, \
                                       convergence_limit = 1.0e-4):
    '''
    estimate_waterdepth_from_discharge :
                       function to estimate the water depth correspondent to a discharge;
                       values are obtained by iteration, starting by a water depth of the
                       previous time-step
    
    input:
    =====
    discharge         : PCRaster map with discharge values (units: m3/s)
    waterdepth        : PCRaster map with water depth to start iteration (units: m)
    mannings_n        : PCRaster map with Manning's coefficient [m^-1/3*s]
    channel_width     : PCRaster map with the width for a rectangular channel (units: m)
    channel_gradient  : PCRaster map with the gradient along the channel (units: m/m)
    beta              : float, exponent in equatio Area = alpha * Q ** beta
    max_iterations    : integer, maximum number of iterations to calculate the correspondent
                        water depth to the input discharge
    convergence_limit : float, minimum convergence value expected to stop the iterative
                        process; calculated as the difference between the currently
                        estimated water depth and the one from the previous time-step
    
    output:
    ======
    waterdepth       : water depth at the end of the time-step based on the discharge
                       and the channel parameters (units: m)
    '''
    
    # set the mask where water is present (units: m3/s)
    discharge_mask = discharge > 0
    
    # set the number of iterations and convergence
    icnt = 0
    convergence = False
    
    # and iterate over all time steps
    while icnt < max_iterations and not convergence:

        # set the old water 
        waterdepth_old = pcr.max(0.001, waterdepth)

        # get the wetted perimeter [m] for a rectangular channel and the 
        # corresponding alpha of the equation A = alpha * Q ** beta
        wetted_perimeter = channel_width + 2 * waterdepth_old
        alpha = (mannings_n * wetted_perimeter ** (2.0/3.0) * channel_gradient ** -0.5) ** beta

        # compute the new water table
        wetted_area = alpha * discharge ** beta
        waterdepth = pcr.max(0.001, wetted_area / channel_width)
        
        # compare the water depth
        conv_value  = pcr.cellvalue(pcr.mapmaximum(pcr.abs(waterdepth - waterdepth_old)),1)[0]
        convergence = conv_value < convergence_limit
        
        # update icnt
        icnt = icnt + 1
    
    # set the water depth
    waterdepth = pcr.ifthenelse(discharge_mask, waterdepth, pcr.scalar(0))
    
    # set the message string
    message_str = 'water depth converged after %d iterations with a maximum deviation of %.3g' % (icnt, conv_value)
    print(message_str)
    
    # return water depth
    return waterdepth


###################
# class definition #
###################

"""
water_management: read in the monthly long-term availability for the discharge \
and the base flow.
"""

class water_management(object):

    """
water management module: class that holds the water management module of \
the QUAlloc model and manages water use through allocation and withdrawal.

Currently, only surface water and groundwater are considered but this can be \
expanded to include desalinated water too.
Also, it assumed that each cell receives water from a predined set of points \
for one resource only, i.e., zones are mutually exclusive per resource.
At the moment, all withdrawal points are fixed in time and the withdrawal \
capacity is not considered.
time_increment                    : time increment to be used;
                                    currently only monthly and yearly allowed
                                    and this refers to the update of the with-
                                    drawals given the long-term availability;
desalwater_allocation_zones,
surfacewater_allocation_zones,
groundwater_allocation_zones      : nominal map with the ID of the allocation
                                    zone per resource identified by non-zero
                                    values;
desalwater_withdrawal_points,
surfacewater_withdrawal_points,
groundwater_withdrawal_points     : map with the withdrawal points falling 
                                    in a zone; this should be an ordinal map
                                    that gives each point a unique ID and
                                    that can be used to track local with-
                                    drawals and assign the withdrawal capac-
                                    ity.
surfacewater_withdrawal_capacity,
groundwater_withdrawal_capacity   : the capacity [m^3/day] to withdraw water
                                    from the resource specified. Currently,
                                    this can only be None, in which case the
                                    abstraction is unlimited or set to a pre-
                                    defined rate at the location. If with-
                                    drawal points were to become dynamic,
                                    this could be linked to the allocated
                                    demand and capped by a pre-defined capac-
                                    ity per region, as currenly is done in
                                    PCR-GLOBWB 2 with the groundwater pumping
                                    capacity.
surfacewater_update_weight,
groundwater_update_weight         : weight, defined as 1/N years, by which
                                    the long-term average water availability,
                                    abstractions and return flows are updat-
                                    ed with the current conditions [year^-1],
                                    weight should be greater than 0 and less
                                    than 1;
prioritization                    : priority for the allocation of the total
                                    demand to the available withdrawals, dimen-
                                    sionless, as a dictionary with the sector
                                    names as keys and scalar PCRaster fields
                                    as values; default value is None, in which
                                    case all sectors are given equal weight;
                                    priorities are given in ascending order,
                                    the lowest value having the highest prior-
                                    ity; priorities can be defined locally as
                                    maps but are spatially compared; currently,
                                    a single prioritization is provided but
                                    internally the sources surface water and
                                    groundwater are already distinguished.

Long-term averages per month need to be provided for the following variables
as monthly netCDFs:
surfacewater_longterm_dischage_ini        : long-term average monthly surface 
                                            water discharge [m3/s]
surfacewater_longterm_runoff_ini           : long-term average monthly surface 
                                            water runoff [m/day]
groundwater_longterm_storage_ini          : long-term average monthly groundwater
                                            storage [m per day]

# and the total return flow needs to be initialized, default is None,
# in which case it is set to zero
total_return_flow_ini                      : total return flow [m3/day]

"""
    
    def __init__(self, \
        landmask, \
        time_increment, \
        cellarea, \
        time_step_length, \
        desalwater_allocation_zones, \
        desalwater_withdrawal_points, \
        groundwater_allocation_zones, \
        groundwater_withdrawal_points, \
        groundwater_withdrawal_capacity, \
        groundwater_update_weight, \
        groundwater_longterm_storage, \
        groundwater_longterm_potential_withdrawal, \
        surfacewater_allocation_zones, \
        surfacewater_withdrawal_points, \
        surfacewater_withdrawal_capacity, \
        surfacewater_update_weight, \
        surfacewater_longterm_discharge, \
        surfacewater_longterm_runoff, \
        surfacewater_longterm_potential_withdrawal, \
        total_return_flow_ini = None, \
        prioritization       = None, \
        sector_names         = ['irrigation','domestic','industry','livestock'], \
        source_names         = ['groundwater','surfacewater'], \
        withdrawal_names     = ['renewable','nonrenewable'], \
        water_quality_flag                 = False, \
        desalinated_water_use_flag         = False, \
        groundwater_pumping_capacity_flag  = False, \
        surfacewater_pumping_capacity_flag = False, \
        ):

        '''
water management class requires the following input for its initialization:
    
    Input:
    ======
        landmask,
        time_increment,
        groundwater_allocation_zones,
        groundwater_withdrawal_points,
        groundwater_withdrawal_capacity,
        groundwater_update_weight,
        groundwater_longterm_baseflow,
        surfacewater_allocation_zones,
        surfacewater_withdrawal_points,
        surfacewater_withdrawal_capacity,
        surfacewater_update_weight,
        surfacewater_longterm_discharge,
        surfacewater_longterm_runoff,
        total_return_flow_ini = None
        prioritization       = None
        sector_names         = ['irrigation','domestic','industry','livestock']
        source_names         = ['groundwater','surfacewater']
        water_quality_flag                 = False
        desalinated_water_use_flag         = False
        groundwater_pumping_capacity_flag  = False
        surfacewater_pumping_capacity_flag = False
See doc string of class for detailed info.
  
'''
        
        # initialize the object
        object.__init__(self)
        
        # set model land mask
        self.landmask = landmask
        
        # source names to be used to process the withdrawals
        # NOTE: names could be used for processing but is largely ignored here
        # as the processing of the different sources may vary.
        self.source_names = source_names
        
        # sector names used to process the demand
        self.sector_names = sector_names
        
        # withdrawal types; renewable and non-renewable; this is hard-coded
        self.withdrawal_names = withdrawal_names
        
        # set the information on the settings to allocate demand to the avail-
        # ability; these are currently set hard coded here and refer to the
        # allocation functions
        # time increment to be used; currently only monthly and yearly allowed
        # and this refers to the update of the withdrawals given the 
        # long-term availability
        self.time_increment   = time_increment
        self.time_step_length = time_step_length
        
        # use_local_first is a boolean PCRaster map that indicates if the local
        # availability should be used first; the flag is a boolean identifying
        # its overall setting to True or False for logging purposes
        self.use_local_first    = pcr.spatial(pcr.boolean(1))
        self.use_local_first_flag = pcr.cellvalue(pcr.mapmaximum(pcr.scalar( \
                                                  self.use_local_first)), 1)[0] == 1
        
        # reallocate_surplus is a boolean variable (True, False) indicating
        # that any surplus will be used to satisfy any local demand.
        self.reallocate_surplus = True
        
        # source_unmet_demand: either any or all of the available sources;
        # formerly, by default set to both surface water and groundwater that
        # are then allocated proportionally to the withdrawals;
        # note that this is a sequence of the source names
        # and that the first is used as default
        self.sources_unmet_demand = ['surfacewater','groundwater']
        
        # pumping capacity flag is a boolean value to determine if groundwater
        # withdrawal capacity is considered to cap the groundwater abstraction
        self.pumping_capacity_flag = \
                          {'surfacewater': surfacewater_pumping_capacity_flag, \
                           'groundwater' : groundwater_pumping_capacity_flag}
        
        # desalinated water use flag is a boolean value to determine if desalinated
        # water use is considered to contribute tothe water supply before the
        # allocation from surface water and groundwater
        self.desalinated_water_use_flag = desalinated_water_use_flag
        
        # surface water quality flag is a boolean value to determine if surface
        # water quality is considered to cap the surface water availability
        self.water_quality_flag = water_quality_flag
        
        # set the information about withdrawal and allocation, this includes:
        # the ID per allocation zone for which the available water 
        # is matched to the demand;
        # the withdrawal points, where water will be extracted, and that are
        # identified by the ID of the allocation zones
        # the withdrawal capacity at the withdrawal points that limit the ab-
        # straction in [m3/day].
        # Currently, all entries are maps (nomimal for the IDs, scalar for the
        # withdrawal capacity). When the withdrawal capacity is set to None,
        # the withdrawal is unlimited.
        self.groundwater_withdrawal_points    = groundwater_withdrawal_points
        self.groundwater_withdrawal_capacity  = groundwater_withdrawal_capacity
        self.groundwater_update_weight        = groundwater_update_weight
        
        self.surfacewater_withdrawal_points   = surfacewater_withdrawal_points
        self.surfacewater_withdrawal_capacity = surfacewater_withdrawal_capacity
        self.surfacewater_update_weight       = surfacewater_update_weight  
        
        self.desalwater_withdrawal_points     = desalwater_withdrawal_points
        
        # mask extension of withdrawal capacity
        if not isinstance(self.groundwater_withdrawal_capacity, NoneType):
            self.groundwater_withdrawal_capacity = \
                 pcr.ifthen(self.landmask, self.groundwater_withdrawal_capacity)
        if not isinstance(self.surfacewater_withdrawal_capacity, NoneType):
            self.surfacewater_withdrawal_capacity = \
                 pcr.ifthen(self.landmask, self.surfacewater_withdrawal_capacity)
        
        # [ allocation zones ]
        # define zones per sector
        self.groundwater_allocation_zones  = dict((sector_name, groundwater_allocation_zones) \
                                                 for sector_name in sector_names)
        self.surfacewater_allocation_zones = dict((sector_name, surfacewater_allocation_zones) \
                                                  for sector_name in sector_names)
        self.desalwater_allocation_zones  =  dict((sector_name, desalwater_allocation_zones) \
                                                 for sector_name in sector_names)
        
        # define the sectors that can use water from desalinated source
        self.sector_names_desalwater = ['domestic','industry','manufacture']
        
        # define the sectors that can withdraw water only from surfacewater and
        # from its local zone (i.e., same cell)
        self.sectors_local_surfacewater = ['thermoelectric','environment']
        
        # update the allocation zone the sectors where water can be withdrawn
        # only from local zone, thus, each pixel is its own allocation zone
        for sector_name in self.sectors_local_surfacewater:
            if sector_name in self.sector_names:
                self.surfacewater_allocation_zones[sector_name] = \
                        pcr.ifthen(pcr.scalar(surfacewater_allocation_zones)>0, \
                                   pcr.nominal(pcr.uniqueid(pcr.ifthen(pcr.scalar(surfacewater_allocation_zones)>0, \
                                                                       pcr.boolean(1)))))
                
                self.groundwater_allocation_zones[sector_name]  = \
                        pcr.ifthen(pcr.scalar(groundwater_allocation_zones)>0, \
                                   pcr.nominal(pcr.uniqueid(pcr.ifthen(pcr.scalar(groundwater_allocation_zones)>0, \
                                                                                  pcr.boolean(1)))))
                
                # alternative, make all groundwater allocation zones equal to missing values
                # but, watch out! all later calculations should be covered with zeros
        
        # [ prioritization ]
        # a dictionary with source and sector; if not specified (None) default is
        # set to an equal weight; otherwise, priorities are given in ascending
        # order, the lowest value having the highest priority; priorities can be 
        # defined locally as maps but are spatially compared
        source_names_prioritization = deepcopy(self.source_names)
        if desalinated_water_use_flag and 'desalwater' not in source_names:
            source_names_prioritization.append('desalwater')
        
        self.prioritization = dict((source_name, \
                                    dict((sector_name, \
                                          pcr.ifthen(self.landmask, \
                                                     pcr.spatial(pcr.scalar(1)))) \
                                         for sector_name in self.sector_names)) \
                                   for source_name in source_names_prioritization)
        
        if isinstance(prioritization, dict):
            for source_name in source_names_prioritization:
                for sector_name in self.sector_names:
                    self.prioritization[source_name][sector_name] = \
                        pcr.cover(prioritization[source_name][sector_name], \
                                  self.prioritization[source_name][sector_name])
        
        # update the prioritization including the effect of the allocation zones;
        # sector that can withdraw water only from local sources has higher priority
        # note: it only applies to surfacewater and groundwater sources
        for source_name in self.source_names:
            zones = getattr(self, '%s_allocation_zones' % source_name)
            for sector_name in self.sector_names:
                n_cells = get_zonal_total( \
                                          pcr.ifthen(pcr.defined(zones[sector_name]), \
                                                     pcr.scalar(1)), \
                                                     zones[sector_name])
                self.prioritization[source_name][sector_name] = \
                     self.prioritization[source_name][sector_name] * n_cells
        
        # [ long-term ]
        # set long-term variables use to calculate the long-term availability
        # groundwater_longterm_storage    (units: m per day)
        # surfacewater_longterm_discharge (units: m3/s)
        # surfacewater_longterm_runoff     (units: m/day)
        self.groundwater_longterm_storage     = groundwater_longterm_storage
        self.surfacewater_longterm_discharge  = surfacewater_longterm_discharge
        self.surfacewater_longterm_runoff      = surfacewater_longterm_runoff
        
        # set long-term variables use to define the long-term potential withdrawals
        # note:
        #     if values are unknown, long-term groundwater storage and long-term surface water
        #     total runoff could be used to initialize the variables (in volume over time)
        # (units: m3/day)
        self.groundwater_longterm_potential_withdrawal  = groundwater_longterm_potential_withdrawal
        self.surfacewater_longterm_potential_withdrawal = surfacewater_longterm_potential_withdrawal
        
        # set a list of sorted dates of long-term water availability
        self.groundwater_longterm_storage_dates  = \
                         sorted(list(self.groundwater_longterm_storage.keys()))
        self.surfacewater_longterm_discharge_dates = \
                         sorted(list(self.surfacewater_longterm_discharge.keys()))
        self.surfacewater_longterm_runoff_dates = \
                         sorted(list(self.surfacewater_longterm_runoff.keys()))
        
        self.groundwater_longterm_pot_withdrawal_dates = \
                         sorted(list(self.groundwater_longterm_potential_withdrawal.keys()))
        self.surfacewater_longterm_pot_withdrawal_dates = \
                         sorted(list(self.surfacewater_longterm_potential_withdrawal.keys()))
        
        # get the total (annual average) long-term availability
        # they keep the same units as their correspondent monthly long-term counterparts
        self.groundwater_total_storage    = pcr.scalar(0)
        self.surfacewater_total_discharge = pcr.scalar(0)
        self.surfacewater_total_runoff     = pcr.scalar(0)
        self.update_total_water_availability()
        
        # [ initialize outputs ]
        # set the sectoral gross and net water demand as volume per cell
        self.gross_demand = dict((sector_name, \
                                  pcr.spatial(pcr.scalar(0))) \
                                 for sector_name in self.sector_names)
        self.net_demand   = dict((sector_name, \
                                  pcr.spatial(pcr.scalar(0))) \
                                 for sector_name in self.sector_names)
        
        # set the total gross and net demand per cell
        # as well as the consumption and return flows,
        # the total withdrawal and total allocated demand
        self.total_gross_demand = pcr.spatial(pcr.scalar(0))
        self.total_net_demand   = pcr.spatial(pcr.scalar(0))
        self.total_consumption  = pcr.spatial(pcr.scalar(0))
        self.total_return_flow   = pcr.spatial(pcr.scalar(0))
        self.total_withdrawal   = pcr.spatial(pcr.scalar(0))
        self.total_allocation   = pcr.spatial(pcr.scalar(0))
        
        # use the input for the total return flow if provided
        # (units: m3/day)
        if not isinstance(total_return_flow_ini, NoneType):
            self.total_return_flow = pcr.cover(pcr.scalar(total_return_flow_ini), \
                                               self.total_return_flow)
        
        # initialize the potetial renewable and non-renewable withdrawal for the sources
        self.potential_renewable_withdrawal = \
                                          dict((source_name, \
                                                pcr.spatial(pcr.scalar(0))) \
                                          for source_name in self.source_names)
        self.potential_nonrenewable_withdrawal = \
                                          dict((source_name, \
                                                pcr.spatial(pcr.scalar(0))) \
                                          for source_name in self.source_names)
        
        # idem for the actual ones
        self.actual_renewable_withdrawal = \
                                          dict((source_name, \
                                                pcr.spatial(pcr.scalar(0))) \
                                          for source_name in self.source_names)
        self.actual_nonrenewable_withdrawal = \
                                          dict((source_name, \
                                                pcr.spatial(pcr.scalar(0))) \
                                          for source_name in self.source_names)
        
        # idem for the unused withdrawal
        self.unused_renewable_withdrawal = \
                                          dict((source_name, \
                                                pcr.spatial(pcr.scalar(0))) \
                                          for source_name in self.source_names)
        self.unused_nonrenewable_withdrawal = \
                                          dict((source_name, \
                                                pcr.spatial(pcr.scalar(0))) \
                                          for source_name in self.source_names)
        
        # allocated withdrawal and demand, organized as nested dictionaries
        # with combined withdrawal / source names as key of a dictionary
        # with the sector names as key, as well as the consumption and 
        # return flow at the cell level on the basis of the allocation per cell
        self.consumed_demand_per_sector      = {}
        self.return_flow_demand_per_sector    = {}
        self.allocated_withdrawal_per_sector = {}
        self.allocated_demand_per_sector     = {}
        
        for withdrawal_name in self.withdrawal_names:
            for source_name in self.source_names:
                
                key = get_key([withdrawal_name, source_name])
                
                self.allocated_withdrawal_per_sector[key] = dict((sector_name, \
                                                                  pcr.spatial(pcr.scalar(0))) \
                                                            for sector_name in sector_names)
                
                self.allocated_demand_per_sector[key]     = dict((sector_name, \
                                                                  pcr.spatial(pcr.scalar(0))) \
                                                            for sector_name in sector_names)
                
                self.consumed_demand_per_sector[key]      = dict((sector_name, \
                                                                  pcr.spatial(pcr.scalar(0))) \
                                                            for sector_name in sector_names)
                
                self.return_flow_demand_per_sector[key]   = dict((sector_name, \
                                                                  pcr.spatial(pcr.scalar(0))) \
                                                             for sector_name in sector_names)
        
        # idem for desalinated water use
        self.allocated_withdrawal_per_sector_desalwater = dict((sector_name, \
                                                                pcr.spatial(pcr.scalar(0))) \
                                                           for sector_name in sector_names)
        
        self.allocated_demand_per_sector_desalwater     = dict((sector_name, \
                                                                pcr.spatial(pcr.scalar(0))) \
                                                          for sector_name in sector_names)
        
        self.consumed_demand_per_sector_desalwater      = dict((sector_name, \
                                                                pcr.spatial(pcr.scalar(0))) \
                                                          for sector_name in sector_names)
        
        self.return_flow_demand_per_sector_desalwater    = dict((sector_name, \
                                                                pcr.spatial(pcr.scalar(0))) \
                                                           for sector_name in sector_names)
        
        # [ states ]
        # set the state names
        self.report_state_info = { \
                  'total_return_flow'                           : 'total_return_flow', \
                  'surfacewater_longterm_discharge'            : 'surfacewater_longterm_discharge', \
                  'surfacewater_longterm_runoff'                : 'surfacewater_longterm_runoff', \
                  'groundwater_longterm_storage'               : 'groundwater_longterm_storage', \
                  'groundwater_longterm_potential_withdrawal'  : 'groundwater_longterm_potential_withdrawal', \
                  'surfacewater_longterm_potential_withdrawal' : 'surfacewater_longterm_potential_withdrawal'}
        
        #for source_name in self.source_names:    # <------------------------------------------------- pumping capacity
        #    if self.pumping_capacity_flag[source_name]:
        #        if source_name == 'groundwater':
        #            self.report_state_info['groundwater_longterm_potential_withdrawal'] = \
        #                                   'groundwater_longterm_potential_withdrawal'
        #        if source_name == 'surfacewater':
        #            self.report_state_info['surfacewater_longterm_potential_withdrawal'] = \
        #                                   'surfacewater_longterm_potential_withdrawal'
        
        # [ reporting ]
        # create message string on the selected processing
        message_str = 'water demand options are processed with the following options:'
        for option_str in ['time_increment', 'use_local_first_flag', 'reallocate_surplus']:
            message_str = str.join('\n\t', \
                                   (message_str, \
                                    '%-20s: %s' % (option_str, \
                                                  getattr(self, option_str))))
        if isinstance(self.sources_unmet_demand, NoneType):
            message_str = str.join('\n\t', \
                                   (message_str, \
                                    'No source is set for non-renewable',\
                                    'withdrawals and any demand that', \
                                    'is not met by the renewable resources,', \
                                    'is forfeited.'))
        else:
            sub_message_str = str.join(', ', \
                                      (self.source_names))
            message_str     = str.join('\n\t', \
                                      (message_str, \
                                      'The following non-renewable sources are set to satisfy any unmet demand:', \
                                      sub_message_str))
        for source_name in self.source_names:
            if self.pumping_capacity_flag[source_name]:
                message_str = str.join('\n\t', \
                                      (message_str, \
                                      'Pumping capacity is considered to limit the %s abstraction.' % source_name))
        if self.water_quality_flag:
            message_str     = str.join('\n\t', \
                                      (message_str, \
                                      'Water quality is considered to determine the actual water availability.'))
        
        # log message_str
        logger.info(message_str)
        
        # returns None
        return None
    
    
    
    def update_total_water_availability(self):
        '''
        update the total water availability over the year.
        '''
        
        # log message
        logger.info('total water availability updated')
        
        # in the following, the weights are weights for the dates as a dictionary
        # with the latter as key
        # update the long-term total groundwater availability (storage)
        # (units: m per day)
        weights = get_weights_from_dates(self.groundwater_longterm_storage_dates)
        self.groundwater_total_storage  = \
                         sum(list(weights[date] * \
                             self.groundwater_longterm_storage[date] \
                             for date in self.groundwater_longterm_storage_dates))
        
        # update the long-term total surface water availability (discharge)
        # (units: m3/s)
        weights = get_weights_from_dates(self.surfacewater_longterm_discharge_dates)
        self.surfacewater_total_discharge = \
                         sum(list(weights[date] * \
                             self.surfacewater_longterm_discharge[date] \
                             for date in self.surfacewater_longterm_discharge_dates))
        
        # update the long-term total surface water availability (runoff)
        # (units: m/day)
        weights = get_weights_from_dates(self.surfacewater_longterm_runoff_dates)
        self.surfacewater_total_runoff  = \
                         sum(list(weights[date] * \
                             self.surfacewater_longterm_runoff[date] \
                             for date in self.surfacewater_longterm_runoff_dates))
        
        #if self.pumping_capacity_flag['groundwater']:
        #    # update the potential long-term total groundwater withdrawals
        #    weights = get_weights_from_dates(self.groundwater_longterm_pot_withdrawal_dates)
        #    self.groundwater_total_potential_withdrawal  = \
        #                     sum(list(weights[date] * \
        #                         self.groundwater_longterm_potential_withdrawal[date] \
        #                         for date in self.groundwater_longterm_pot_withdrawal_dates))
        #
        #if self.pumping_capacity_flag['surfacewater']:
        #    # update the potential long-term total surface water withdrawals
        #    weights = get_weights_from_dates(self.surfacewater_longterm_pot_withdrawal_dates)
        #    self.surfacewater_total_potential_withdrawal  = \
        #                     sum(list(weights[date] * \
        #                         self.surfacewater_longterm_potential_withdrawal[date] \
        #                         for date in self.surfacewater_longterm_pot_withdrawal_dates))
        
        # returns None
        return None
    
    
    
    def update_withdrawal_capacity(self, \
                                    source_name, \
                                    regional_pumping_limit, \
                                    region_ids, \
                                    region_ratios, \
                                    time_step_length, \
                                    date):
        '''
        update_withdrawal_capacity : 
                                 function to calculate the withdrawal capacity
                                 based on the regional pumping capacity volumes reported
                                 over predefined regions.
        
        input:
        =====
        regional_pumping_limit : PCRaster map with regional values (same value across 
                                 the region) of maximum water pumping capacity
                                 (units: billion cubic meters/year)
        region_ids             : PCRaster map with regional ID (units: nominal)
        region_ratio           : PCRaster map with scalar values that scale the actual
                                 pumping capacity comprised within the land mask (units: m3/m3)
        time_step_length       : integer, number of days in a period (e.g., 30 days/month)
        date                   : string, date of the update
        '''
        
        # scale the regional pumping capacity to the extension of the land mask
        # and convert units of regional pumping limit:
        # from billion cubic meters per year to cubic meters per day
        # (units: m3/year)
        regional_pumping_limit = regional_pumping_limit * region_ratios * 1000000000
        
        # define the long-term potential withdrawal per source
        # (units: m3/day)
        if source_name == 'groundwater':
            longterm_potential_withdrawal = deepcopy(self.groundwater_longterm_potential_withdrawal)
        if source_name == 'surfacewater':
            longterm_potential_withdrawal = deepcopy(self.surfacewater_longterm_potential_withdrawal)
        
        # calculate the rate of water available during the specified month
        # (units: m3/m3)
        maximum_monthly_longterm_potential_withdrawal = \
                                            max_dicts(longterm_potential_withdrawal)
        
        total_maximum_monthly_longterm_potential_withdrawal = \
            get_zonal_total(maximum_monthly_longterm_potential_withdrawal, \
                            region_ids)
        
        withdrawal_rate = \
            pcr_return_val_div_zero(maximum_monthly_longterm_potential_withdrawal, \
                                    total_maximum_monthly_longterm_potential_withdrawal, \
                                    very_small_number)
        
        # update the variable corresponding to the water withdrawal capacity
        # (units: m3/day)
        withdrawal_capacity = regional_pumping_limit * withdrawal_rate / \
                                                      (12 * time_step_length)
        if source_name == 'groundwater':
            self.groundwater_withdrawal_capacity  = withdrawal_capacity
        if source_name == 'surfacewater':
            self.surfacewater_withdrawal_capacity = withdrawal_capacity
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        if verbose:
            dt = f'{str(date.year)[2:]}-{str(date.month).zfill(2)}'
            pcr.report(self.groundwater_withdrawal_capacity, f'{path}/{dt}_groundwater_withdrawal_capacity.map')
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # return None
        return None
    
    
    
    def update_water_demand_for_date(self, \
                                     gross_demand, \
                                     net_demand, \
                                     date):
        
        '''
        update_water_demand_for_date: 
                        function that updates the availability per
                        zone as a function of the date.
        
        input:
        ======
        gross_demand,
        net_demand    : gross and net demand as a dictionary
                        with the sector names as keys (m3/day)
                        as a PCRaster map.
        date          : date of the update.
        '''
        
        # set the sectoral gross and net water demand
        self.gross_demand = dict((sector_name, 
                                  pcr.spatial(pcr.scalar(0))) \
                                 for sector_name in self.sector_names)
        self.net_demand   = dict((sector_name, 
                                  pcr.spatial(pcr.scalar(0))) \
                                 for sector_name in self.sector_names)
        
        # iterate over the sector names and update the values if applicable
        for sector_name in self.sector_names:
            
            # [ gross demand ]
            if sector_name in gross_demand.keys():
                
                # log message
                logger.debug('%s gross water demand set for %s' % \
                             (sector_name, date))
                
                # set the value (units: m3/day)
                self.gross_demand[sector_name] = gross_demand[sector_name]
                
                # add the gross demand to the long-term gross demand for the 
                # present date
                pass
            
            # [ net demand ]
            if sector_name in net_demand.keys():
                
                # log message
                logger.debug('%s net water demand set for %s' % \
                             (sector_name, date))
                
                # set the value(units: m3/day)
                self.net_demand[sector_name] = net_demand[sector_name]
                
                # add the net demand to the long-term net demand for the 
                # present date
                pass
        
        # get the totals: gross and net
        # (units: m3/day)
        self.total_gross_demand = sum_list(list(self.gross_demand.values()))
        self.total_net_demand   = sum_list(list(self.net_demand.values()))
        
        # create an updateable gross demand variable
        # (units: m3/day)
        self.gross_demand_remaining = deepcopy(self.gross_demand)
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        if verbose:
            dt = f'{str(date.year)[2:]}-{str(date.month).zfill(2)}'
            for sector_name in self.sector_names:
                pcr.report(self.gross_demand_remaining[sector_name], f'{path}/{dt}_gross_demand_{sector_name}.map')
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # add the total gross and net demand to the long-term net demand for the 
        # present date
        pass
        
        # log message
        logger.debug('total gross and net demand set for %s' % date)
       
        # returns None
        return None
    
    
    
    def get_longterm_availability_for_date(self, \
                                            date, \
                                            cellarea, \
                                            ldd, \
                                            waterdepth, \
                                            mannings_n, \
                                            channel_gradient, \
                                            channel_width, \
                                            channel_length, \
                                            time_step_seconds = 86400):
        '''
        get_longterm_availability_for_date: 
                                    function that gets the water availability for the time base
                                    to be used (monthly, in which case the date is matched to the 
                                    nearest date in the availability; or yearly, in which case the  
                                    total long-term availability is used)
        
        input:
        =====
        date                      : string, date of the update
        cellarea                  : PCRaster map with area of cell (units: m2)
        ldd                       : PCRaster map with flow directions
        waterdepth                : PCRaster map with water depth to start iteration (units: m)
        mannings_n                : PCRaster map with Manning's coefficient (units: m^-1/3*s)
        channel_gradient          : PCRaster map with the gradient along the channel (units: m/m)
        channel_width             : PCRaster map with the width for a rectangular channel (units: m)
        channel_length            : PCRaster map with the length of the channel (units: m)
        
        output:
        ======
        surfacewater_availability,
        groundwater_availability  : PCRaster map with the long-term surfacewater nad groundwater
                                    availability (m3/day)
        '''
        
        # set the message string to log the information
        message_str = 'Long-term water availability for %s base on %s time increment' % \
                      (date, self.time_increment)
        
        # [ get water availability ] ........................................................................................
        # 1st set values from long-term availability
        # 2nd patch NaNs with pumping capacity (if available)
        # 3rd fill NaNs with zeros
        # 4th limit the availability to pumping capacity (if available) and withdrawal points
        
        # set long-term water availability
        # and panic if it is a missing value ;-p
        if self.time_increment == 'monthly':
            # get the groundwater and surface water availabilty for the matching date
            
            # [ groundwater ]
            # get the time step to update the groundwater availability:
            # groundwater_storage (units: m per day)
            date_index, matched_date, sub_message_str = \
                    match_date_by_julian_number(date, \
                                                self.groundwater_longterm_storage_dates)
            groundwater_storage = self.groundwater_longterm_storage[matched_date]
            
            # add the message on the matching date to the string
            message_str = str.join('\n', \
                                   (message_str, sub_message_str))
            
            # [ surface water ]
            # get the time step to update the surface water availability:
            # surfacewater_discharge (units: m3/s)
            # surfacewater_runoff     (units: m/day)
            date_index, matched_date, sub_message_str = \
                    match_date_by_julian_number(date, \
                                                self.surfacewater_longterm_discharge_dates)
            surfacewater_discharge = self.surfacewater_longterm_discharge[matched_date]
            
            date_index, matched_date, sub_message_str = \
                    match_date_by_julian_number(date, \
                                                self.surfacewater_longterm_runoff_dates)
            surfacewater_runoff = self.surfacewater_longterm_runoff[matched_date]
            
            # add the message on the matching date to the string
            message_str = str.join('\n', \
                                   (message_str, sub_message_str))
        
        elif self.time_increment == 'yearly':
            # set the long-term total availability
            # groundwater_storage    (units: m per day)
            # surfacewater_discharge (units: m3/s)
            # surfacewater_runoff     (units: m/day)
            groundwater_storage    = self.groundwater_total_storage
            surfacewater_discharge = self.surfacewater_total_discharge
            surfacewater_runoff     = self.surfacewater_total_runoff
        
        else:
            logger.error('the option %s for the time increment in the water management module is not allowed!')
            sys.exit()
        
        # get the surface and groundwater availability for date
        #
        # [ surface water ]
        # get the average daily discharge by adding the discharge from
        # the cell upstream and the total runoff of the same cell
        # (units: m3/s)
        discharge = pcr.upstream(ldd, surfacewater_discharge) + \
                    surfacewater_runoff * cellarea / time_step_seconds
        
        # get the average daily water depth corresponding to this discharge
        # (units: m per day)
        channel_depth = estimate_waterdepth_from_discharge( \
                                       discharge        = discharge, \
                                       waterdepth       = waterdepth, \
                                       mannings_n       = mannings_n, \
                                       channel_width    = channel_width, \
                                       channel_gradient = channel_gradient)
        
        # get the surface water availability (units: m3/day)
        surfacewater_availability = channel_depth * channel_width * channel_length
        
        # [ groundwater ]
        # get the groundwater availability (units: m3/day)
        groundwater_availability = groundwater_storage * cellarea
        
        # patch to exclude non-zero availability
        groundwater_availability  = pcr.ifthen(groundwater_availability  >= 0, \
                                               groundwater_availability)
        surfacewater_availability = pcr.ifthen(surfacewater_availability >= 0, \
                                               surfacewater_availability)
        
        # update the availability (depending on the pumping capacity)
        # check on possible missing data; note that updates on the local variable
        # water availability do not interact with the long-term availability
        # stored internally to the instance of the water management module
        missing_availability = self.landmask & (pcr.pcrnot(pcr.defined(groundwater_availability)) |
                                           pcr.pcrnot(pcr.defined(surfacewater_availability)))
        missing_availability_flag = pcr.cellvalue(pcr.mapmaximum(pcr.scalar( \
                                                  missing_availability)), 1)[0] == 1
        
        # if cells are missing, first patch with the pumping capacity if available
        # [ groundwater ]
        if not isinstance(self.groundwater_withdrawal_capacity, NoneType):
            groundwater_availability = pcr.cover(groundwater_availability, \
                                                 self.groundwater_withdrawal_capacity)
        # [ surface water ]
        if not isinstance(self.surfacewater_withdrawal_capacity, NoneType):
            surfacewater_availability = pcr.cover(surfacewater_availability, \
                                                  self.surfacewater_withdrawal_capacity)
        
        # set the warning on the flag of missing availability and if no capacity is defined
        if missing_availability_flag:
            
            # add the warning to the message string
            message_str = str.join('\n', \
                                   (message_str, \
                                    'Warning: long-term water availability contains missing values;'))
            
            if not isinstance(self.groundwater_withdrawal_capacity, NoneType):
                message_str = str.join('\n', \
                                       (message_str, \
                                        'by default it is set to the pumping capacity for groundwater;'))
            else:
                message_str = str.join('\n', \
                                       (message_str, \
                                        'by default, zero values are added for groundwater;'))
                
            if not isinstance(self.surfacewater_withdrawal_capacity, NoneType):
                message_str = str.join('\n', \
                                       (message_str, \
                                        'by default it is set to the pumping capacity for surface water.'))
            else:
                message_str = str.join('\n', \
                                       (message_str, \
                                        'by default, zero values are added for surface water.'))
        
        # cover the availability with zeros over the land mask
        surfacewater_availability = pcr.ifthen(self.landmask, \
                                               pcr.cover(surfacewater_availability, 0))
        groundwater_availability  = pcr.ifthen(self.landmask, \
                                               pcr.cover(groundwater_availability, 0))
        
        # and limit the availability to the pumping capacity
        if not isinstance(self.groundwater_withdrawal_capacity, NoneType):
            message_str = str.join('\n', \
                                   (message_str, \
                                    'groundwater withdrawals are limited to the pumping capacity'))
            groundwater_availability = pcr.min(self.groundwater_withdrawal_capacity, \
                                                groundwater_availability)
        
        if not isinstance(self.surfacewater_withdrawal_capacity, NoneType):
            message_str = str.join('\n', \
                                   (message_str, \
                                    'surface water withdrawals are limited to the pumping capacity'))
            surfacewater_availability = pcr.min(self.surfacewater_withdrawal_capacity, \
                                                surfacewater_availability)
        
        # and the availability to the withdrawal points
        groundwater_availability  = pcr.ifthenelse(pcr.scalar(self.groundwater_withdrawal_points)  != 0, \
                                                   groundwater_availability,  0)
        surfacewater_availability = pcr.ifthenelse(pcr.scalar(self.surfacewater_withdrawal_points) != 0, \
                                                   surfacewater_availability, 0)
        
        # return the surface and groundwater long-term availability
        return surfacewater_availability, groundwater_availability
    
    
    
    def update_environmental_flow_requirements_for_date(self, \
                                                        surfacewater_availability, \
                                                        surfacewater_depth, \
                                                        mannings_n, \
                                                        channel_width, \
                                                        channel_gradient, \
                                                        channel_length, \
                                                        time_step_seconds = 86400):
        '''
        update_environmental_flow_requirements:
                                    function that updates the environmental flow requirements based on
                                    the actual water depth
        input:
        =====
        surfacewater_availability : PCRaster map with surface water long-term availability as the
                                    channel storage (units: m3/day)
        surfacewater_depth        : PCRaster map with water depth at the start of the time-step (units: m)
        mannings_n                : PCRaster map with Manning's coefficient [m^-1/3*s]
        channel_gradient          : PCRaster map with the gradient along the channel (units: m/m)
        channel_width             : PCRaster map with the width for a rectangular channel (units: m)
        channel_length            : PCRaster map with the length of the channel (units: m)
        time_step_seconds         : integer, number of second in a day (i.e., 86400 sec/d)
        
        output:
        ======
        prioritization            : PCRaster map with prioritization per source and sector and with updated
                                    value for environmental flow demands
        '''
        
        # [ environmental flow requirements ]
        # convert environmental flow requirements units from days to seconds
        # (units: m3/s)
        discharge_environment = self.gross_demand['environment'] / time_step_seconds
        
        # get the water depth correspondent to the environmental flow requirements (units: m)
        channel_depth_environment = estimate_waterdepth_from_discharge( \
                                       discharge        = discharge_environment, \
                                       waterdepth       = surfacewater_depth, \
                                       mannings_n       = mannings_n, \
                                       channel_width    = channel_width, \
                                       channel_gradient = channel_gradient)
        
        # get the volume of environmental flow to be storaged
        # during the time-step and update the variable (units: m3/day)
        self.gross_demand['environment'] = channel_depth_environment * channel_width * channel_length
        self.net_demand['environment'] = deepcopy(self.gross_demand['environment'])
        
        # [ surface water long-term availability ]
        # get the water depth correspondent to the surface water long-term availability 
        # as channel storage (units: m)
        channel_depth = surfacewater_availability / (channel_width * channel_length)
        
        # update the environmental flow prioritization based on the possibility
        # of supplying this demand
        prioritization = deepcopy(self.prioritization)
        prioritization['surfacewater']['environment'] = \
                            self.prioritization['surfacewater']['environment'] * \
                            pcr.max(0.1, \
                                    pcr_return_val_div_zero(channel_depth, \
                                                            channel_depth_environment, \
                                                            very_small_number))
        
        # return prioritization
        return prioritization
    
    
    
    def allocate_desalinated_water_for_date(self, \
                                             availability,
                                             date):
        '''
        desalinated_water_allocation_for_date: 
                       function that updates the potential withdrawal as a function 
                       of the date to extract the water availability and the internal
                       model settings for the time base to be used (monthly, in which
                       case the date is matched to the nearest date in the availability;
                       or yearly, in which case the long-term availability is used) and
                       the allocation settings.
                       Sets the desalinated water use internally and updates the sectoral
                       gross demands.
        
        input:
        =====
        availability : PCRaster maps with desalinated water availability (m3/day)
        date         : date of the update
        '''
        
        # set the message string to log the information
        message_str = 'Desalinated water use for %s.' % (date)
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        if verbose:
            dt = f'{str(date.year)[2:]}-{str(date.month).zfill(2)}'
            pcr.report(availability, f'{path}/{dt}_shortterm_desalwater_availability.map')
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # [ set up input data ] .............................................................................................
        # define starting conditions
        remaining_availability      = deepcopy(availability)
        
        unmet_demand_per_sector     = deepcopy(self.gross_demand)
        sector_names_no_desalwater  = [sector_name for sector_name in self.sector_names \
                                       if sector_name not in self.sector_names_desalwater]
        for sector_name in sector_names_no_desalwater:
            unmet_demand_per_sector[sector_name] = pcr.spatial(pcr.scalar(0))
        
        met_demand_per_sector       = dict((sector_name, \
                                            pcr.ifthen(self.gross_demand[sector_name] >= 0.0, \
                                                       pcr.scalar(0.0))) \
                                           for sector_name in self.sector_names)
        
        withdrawal_per_sector       = dict((sector_name, \
                                            pcr.spatial(pcr.scalar(0.0))) \
                                           for sector_name in self.sector_names)
        
        allocated_demand_per_sector = dict((sector_name, \
                                            pcr.spatial(pcr.scalar(0.0))) \
                                           for sector_name in self.sector_names)
        
        # define allocation zones
        zones_per_sector = self.desalwater_allocation_zones
        zones = self.desalwater_allocation_zones['domestic']
        
        # define suitability masks of water quality per sector
        # considering a perfect water suitability after desalination
        # (suitability = 1)
        suitability_per_sector = dict((sector_name, \
                                       pcr.spatial(pcr.scalar(1))) \
                                      for sector_name in self.sector_names)
        
        # define initial exit conditions
        iter_allocation = 1
        exit_condition  = False
        max_iter_allocation = 25
        
        # [ start water distribution ] ......................................................................................
        # iterate until either the demand is met or the supply is exhausted
        while not exit_condition:
            
            # define initial values of zonal demand and supply
            # (units: m3/day)
            if iter_allocation == 1:
                totz_demand_old = get_zonal_total( \
                                      local_values = pcr.max(0, \
                                                             sum_list(list(unmet_demand_per_sector.values()))), \
                                      zones        = zones)
                
                totz_supply_old = get_zonal_total( \
                                      local_values = remaining_availability, \
                                      zones        = zones)
            else:
                totz_demand_old = deepcopy(total_zonal_demand)
                totz_supply_old = deepcopy(total_zonal_supply)
            
            # define water availability weights per sector accounting for water quality
            source_name = 'desalwater'
            weights_per_sector = \
                self.water_quality.get_weights_availability_per_sector( \
                         source_name               = source_name, \
                         sector_names              = self.sector_names, \
                         prioritization_per_sector = self.prioritization[source_name], \
                         suitability_per_sector    = suitability_per_sector, \
                         demand_per_sector         = unmet_demand_per_sector, \
                         availability              = remaining_availability, \
                         zones_per_sector          = zones_per_sector, \
                         )
            
            # evaluate potential water withdrawal
            # allocate the sectoral gross demand for current day to the demand,
            # return the withdrawal, allocated demand, met and unmet demands per sector
            for sector_name in self.sector_names:
                
                # calculate the assigned available water per sector
                remaining_availability_sector = \
                    remaining_availability * weights_per_sector[sector_name]
                
                # calculate water withdrawal and allocation
                tmp_withdrawal, tmp_allocated_demand, tmp_met_demand, tmp_unmet_demand, sub_message_str = \
                    allocate_demand_to_availability_with_options( \
                        demand             = unmet_demand_per_sector[sector_name], \
                        availability       = {'desalwater' : remaining_availability_sector}, \
                        zones              = {'desalwater' : zones_per_sector[sector_name]}, \
                        source_names       = ['desalwater'], \
                        use_local_first     = self.use_local_first, \
                        reallocate_surplus = self.reallocate_surplus)
                
                # update water withdrawal and demand values
                met_demand_per_sector[sector_name]   = \
                            met_demand_per_sector[sector_name] + tmp_met_demand
                
                unmet_demand_per_sector[sector_name] = \
                    pcr.max(0,
                            unmet_demand_per_sector[sector_name] - tmp_met_demand)
                
                withdrawal_per_sector[sector_name] = \
                    withdrawal_per_sector[sector_name] + tmp_withdrawal['desalwater']
                
                allocated_demand_per_sector[sector_name] = \
                    allocated_demand_per_sector[sector_name] + tmp_allocated_demand['desalwater']
            
            # aggregate withdrawals and allocations by source
            withdrawal       = sum_list(list(withdrawal_per_sector.values()))
            #allocated_demand = sum_list(list(allocated_demand_per_sector.values()))
            
            # update remaining water available
            remaining_availability = pcr.max(0, \
                                             availability - withdrawal)
            
            # update iteration number and exit condition
            total_zonal_demand = get_zonal_total( \
                                     local_values = sum_list(list(unmet_demand_per_sector.values())), \
                                     zones        = zones)
            total_zonal_supply = get_zonal_total( \
                                     local_values = remaining_availability, \
                                     zones        = zones)
            
            update_mask     = (total_zonal_demand < totz_demand_old) & \
                              (total_zonal_supply < totz_supply_old)
            
            iter_allocation = iter_allocation + 1
            exit_condition  = (pcr.cellvalue(pcr.mapmaximum(pcr.scalar(update_mask)), 1)[0] == 0) | \
                              (iter_allocation > max_iter_allocation)
        
        # [ ends water distribution ] .......................................................................................
        
        # set potential renewable withdrawals per sector, met demands per sector
        # (units: m3/day)
        self.allocated_withdrawal_per_sector_desalwater = \
                                 dict((sector_name, \
                                       withdrawal_per_sector[sector_name]) \
                                      for sector_name in self.sector_names)
        self.allocated_demand_per_sector_desalwater = \
                                 dict((sector_name, \
                                       allocated_demand_per_sector[sector_name]) \
                                      for sector_name in self.sector_names)
        
        # update gross sectoral demands substracting the met demand per sector
        # using desalinated water;
        # Note: met demands for irrigation, thermoelectric and environment are zero
        # (units: m3/day)
        for sector_name in self.sector_names:
            self.gross_demand_remaining[sector_name] = \
                  pcr.ifthenelse(self.gross_demand[sector_name] - met_demand_per_sector[sector_name] > 0, \
                                 self.gross_demand[sector_name] - met_demand_per_sector[sector_name], \
                                 pcr.scalar(0))
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        if verbose:
            pcr.report(remaining_availability, f'{path}/{dt}_shortterm_desalwater_availability_remaining.map')
            for sector_name in self.sector_names:
                pcr.report(allocated_demand_per_sector[sector_name], f'{path}/{dt}_allocated_demand_desalwater_{sector_name}.map')
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # update the message str
        message_str = str.join('\n', \
                               (message_str, sub_message_str))
        
        # return None
        return None
    
    
    
    
    def update_potential_withdrawals_for_date(self, \
                                               availability, \
                                               prioritization, \
                                               date):
        '''
        update_potential_withdrawals_for_date: 
                       function that updates the potential withdrawal as a function 
                       of the date to extract the water availability and the internal
                       model settings for the time base to be used (monthly, in which
                       case the date is matched to the nearest date in the availability;
                       or yearly, in which case the long-term availability is used) and
                       the allocation settings.
                       Sets the potential withdrawal internally and returns it; the
                       similar approach is followed for the actual withdrawals.
        
        input:
        =====
        availability : dictionary with source names (string) as keys and PCRaster maps with
                       the long-term surfacewater and groundwater availability (m3/day)
        date         : date of the update
        
        output:
        ======
        potential_renewable_withdrawal,
        potential_nonrenewable_withdrawal:
                       dictionary of the potential withdrawal as volume over the period 
                       per cell per source
        '''
        
        # set the message string to log the information
        message_str = 'Potential water withdrawals estimated for %s based on the %s long-term availability.' % \
                      (date, self.time_increment)
        
        # [ initialize variables ] ..........................................................................................
        # potetial renewable and non-renewable withdrawal for the sources
        self.potential_renewable_withdrawal = \
                                          dict((source_name, \
                                                pcr.ifthen(self.landmask, \
                                                           pcr.spatial(pcr.scalar(0)))) \
                                          for source_name in self.source_names)
        self.potential_nonrenewable_withdrawal = \
                                          dict((source_name, \
                                                pcr.ifthen(self.landmask, \
                                                           pcr.spatial(pcr.scalar(0)))) \
                                          for source_name in self.source_names)
        self.potential_renewable_withdrawal_per_sector = \
                                          dict((source_name, \
                                                dict((sector_name, \
                                                      pcr.ifthen(self.landmask, \
                                                                 pcr.spatial(pcr.scalar(0)))) \
                                                     for sector_name in self.sector_names)) \
                                               for source_name in self.source_names)
        self.potential_nonrenewable_withdrawal_per_sector = \
                                          dict((source_name, \
                                                dict((sector_name, \
                                                      pcr.ifthen(self.landmask, \
                                                                 pcr.spatial(pcr.scalar(0)))) \
                                                     for sector_name in self.sector_names)) \
                                               for source_name in self.source_names)
        
        # actual renewable and non-renewable withdrawal for the sources
        self.actual_renewable_withdrawal = \
                                          dict((source_name, \
                                                pcr.spatial(pcr.scalar(0))) \
                                          for source_name in self.source_names)
        self.actual_nonrenewable_withdrawal = \
                                          dict((source_name, \
                                                pcr.spatial(pcr.scalar(0))) \
                                          for source_name in self.source_names)
        self.actual_renewable_withdrawal_per_sector = \
                                          dict((source_name, \
                                                dict((sector_name, \
                                                      pcr.spatial(pcr.scalar(0))) \
                                                     for sector_name in self.sector_names)) \
                                               for source_name in self.source_names)
        self.actual_nonrenewable_withdrawal_per_sector = \
                                          dict((source_name, \
                                                dict((sector_name, \
                                                      pcr.spatial(pcr.scalar(0))) \
                                                     for sector_name in self.sector_names)) \
                                               for source_name in self.source_names)
        
        # [ dictionaries ]
        # create a dictionaries of the allocation zones, withdrawal capacities and
        # withdrawal points
        withdrawal_capacity = {'surfacewater' : self.surfacewater_withdrawal_capacity, \
                               'groundwater'  : self.groundwater_withdrawal_capacity}
        
        withdrawal_points   = {'surfacewater' : self.surfacewater_withdrawal_points, \
                               'groundwater'  : self.groundwater_withdrawal_points}
        
        zones_per_sector    = {'surfacewater' : self.surfacewater_allocation_zones, \
                               'groundwater'  : self.groundwater_allocation_zones}
        
        # [ potential renewable withdrawal ] ................................................................................
        #
        # allocate the total gross demand for the current period to the demand
        # and return the withdrawal, allocated demand and the met and unmet demands
        # (units: m3/day)
        withdrawal_per_sector, met_demand_per_sector, sub_message_str = \
               self.allocate_demand_to_renewable_sources( \
                        availability       = availability, \
                        demand_per_sector  = self.gross_demand_remaining, \
                        zones_per_sector   = zones_per_sector, \
                        prioritization     = prioritization, \
                        use_local_first     = self.use_local_first, \
                        reallocate_surplus = self.reallocate_surplus, \
                        date               = date)
        
        # set potential renewable withdrawals per sector, met demands per sector
        self.potential_renewable_withdrawal_per_sector = \
                 dict((source_name, \
                       dict((sector_name, \
                             withdrawal_per_sector[source_name][sector_name]) \
                            for sector_name in self.sector_names)) \
                      for source_name in self.source_names)
        
        met_demand_per_sector = dict((sector_name, \
                                      met_demand_per_sector[sector_name]) \
                                     for sector_name in self.sector_names)
        
        unmet_demand_per_sector = dict((sector_name, \
                                        pcr.max(0.0, \
                                                self.gross_demand_remaining[sector_name] - met_demand_per_sector[sector_name])) \
                                       for sector_name in self.sector_names)
        
        # aggregate and set potential renewable withdrawals, met and unmet demands
        met_demand   = sum_list(list(met_demand_per_sector.values()))
        unmet_demand = sum_list(list(unmet_demand_per_sector.values()))
        self.potential_renewable_withdrawal = \
                            dict((source_name, \
                                  sum_list(list(self.potential_renewable_withdrawal_per_sector[source_name].values()))) \
                                 for source_name in self.source_names)
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        if verbose:
            dt = f'{str(date.year)[2:]}-{str(date.month).zfill(2)}'
            pcr.report(availability['surfacewater'], f'{path}/{dt}_longterm_surfacewater_availability_[noquality].map')
            pcr.report(availability['groundwater'],  f'{path}/{dt}_longterm_groundwater_availability_[noquality].map')
            for sector_name in self.sector_names:
                pcr.report(met_demand_per_sector[sector_name],       f'{path}/{dt}_potential_renewable_met_demand_{sector_name}.map')
                pcr.report(unmet_demand_per_sector[sector_name],     f'{path}/{dt}_potential_renewable_unmet_demand_{sector_name}.map')
            
            for source_name in self.source_names:
                pcr.report(self.potential_renewable_withdrawal[source_name], f'{path}/{dt}_potential_renewable_{source_name}.map')
                for sector_name in self.sector_names:
                    pcr.report(self.potential_renewable_withdrawal_per_sector[source_name][sector_name], f'{path}/{dt}_potential_renewable_{source_name}_{sector_name}.map')
            pcr.report(self.gross_demand_remaining['domestic']-met_demand_per_sector['domestic'],f'{path}/{dt}_potential_renewable_unmet_demand_DOMESTIC.map')
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # update the message str
        message_str = str.join('\n', \
                               (message_str, sub_message_str))
        
        # [ potential non-renewable withdrawal ] ............................................................................
        #
        # update the potential non-renewable withdrawal with the unmet
        # demand from groundwater only according to the potential
        # sources set in sources_unmet_demand for non-renewable withdrawals
        # if no source for unsustainable abstraction is chosen,
        # only sustainable withdrawals are allowed
        # (units: m3/day)
        source_name = 'groundwater'
        self.allocate_unmet_demand_to_nonrenewable_sources( \
                                    unmet_demand_per_sector = unmet_demand_per_sector,\
                                    zones_per_sector        = zones_per_sector, \
                                    withdrawal_capacity     = withdrawal_capacity, \
                                    withdrawal_points       = withdrawal_points)
        
        # [ DELETEME ] verbose <------------------------------------------------------------------------------------------------------------------
        if verbose:
            for source_name in self.source_names:
                pcr.report(self.potential_nonrenewable_withdrawal[source_name], f'{path}/{dt}_potential_nonrenewable_{source_name}.map')
                for sector_name in self.sector_names:
                    pcr.report(self.potential_nonrenewable_withdrawal_per_sector[source_name][sector_name], f'{path}/{dt}_potential_nonrenewable_{source_name}_{sector_name}.map')
        # ----------------------------------------------------------------------------------------------------------------------------------------
        
        # add the information to the message str
        sub_message_str = str.join(' ', \
                                   ('unmet demand is allocated over the available resources', \
                                    'and an average unmet demand of %g [m3] is met by', \
                                    'a withdrawal from non-renewable resources of', \
                                    '%g [m3].'))
          
        sub_message_str = sub_message_str % \
                          (pcr_get_statistics(unmet_demand)['average'], \
                            pcr_get_statistics(sum_list(list( \
                                              self.potential_nonrenewable_withdrawal.values())))['average'])
        
        message_str = str.join('\n', \
                               (message_str, sub_message_str))
        
        # [ potential withdrawal ] ..........................................................................................
        #
        # set the ideal surface water and groundwater potential withdrawal
        # later used to update the long-term potential withdrawal
        # note:
        #      unmet_demand_per_sector['thermoelectric'] was popped out in
        #      allocated_unmet_demand_to_nonrenewable_sources
        # (units: m3/day)
        self.surfacewater_potential_estimated_withdrawal = \
                             deepcopy(self.potential_renewable_withdrawal['surfacewater'])
        self.groundwater_potential_estimated_withdrawal  = \
                             self.potential_renewable_withdrawal['groundwater'] + \
                             sum_list(list(unmet_demand_per_sector.values()))
        
        
        
        #if date.year == 2018 and date.month == 3:
        #    dt = f'{str(date.year)[2:]}-{str(date.month).zfill(2)}'
        #    pcr.report(self.total_gross_demand, f'{path}/{dt}_total_gross_demand.map')
        #    pcr.report(self.potential_renewable_withdrawal['groundwater'], f'{path}/{dt}_potential_renewable_withdrawal_groundwater.map')
        #    pcr.report(self.potential_renewable_withdrawal['surfacewater'], f'{path}/{dt}_potential_renewable_withdrawal_surfacewater.map')
        #    pcr.report(availability['surfacewater'], f'{path}/{dt}_longterm_available_surfacewater.map')
        #    pcr.report(availability['groundwater'],  f'{path}/{dt}_longterm_available_groundwater.map')
        #    pcr.report(sum_list(list(unmet_demand_per_sector.values())), f'{path}/{dt}_unmet_demand.map')
        #    pcr.report(self.groundwater_potential_estimated_withdrawal, f'{path}/{dt}_groundwater_potential_estimated_withdrawal.map')
        #    
        #    #pcr.aguila(self.total_gross_demand, \
        #    #           self.potential_renewable_withdrawal['groundwater'], \
        #    #           sum_list(list(unmet_demand_per_sector.values())), \
        #    #           self.potential_renewable_withdrawal['surfacewater'], \
        #    #           self.groundwater_potential_estimated_withdrawal, \
        #    #           )
        #    pietje
        
        
        # log the message
        logger.info(message_str)
        
        # returns none
        return None
    
    
    
    def allocate_demand_to_renewable_sources(self, \
                                              availability, \
                                              demand_per_sector, \
                                              zones_per_sector, \
                                              prioritization, \
                                              use_local_first, \
                                              reallocate_surplus, \
                                              date):
        '''
        allocate_demand_to_renewable_sources:
                                function to calculate the potential renewable water withdrawals
                                per source and sector for a certain date considering the water quality
        
        input:
        =====
        availability          : dictionary with source names (string) as keys and PCRaster maps with
                                long-term water availability per source as values (units: m3/period)
        demand_per_sector     : dictionary with sector names (string) as keys and PCRaster maps with
                                sectoral gross water demands as values (units: m3/period)
        zones_per_sector      : dictionary with sector names (string) as keys and PCRaster maps with
                                allocation zones per sector (nominal) as values
        prioritization        : dictionary with source names (string) as keys and another dictionary
                                with sector names (string) as keys and PCRaster maps with water use
                                priority in case of sectoral competition as values
        use_local_first        : PCRaster map with boolean values to indicate water is first withdrawn
                                from local source
        reallocate_surplus    : boolean
        
        output:
        ======
        withdrawal_per_sector : dictionary with sector names (string) as keys and PCRaster maps with
                                sum of water withdrawals from renewable surface and groundwater sources
                                per sector as values (units: m3/d)
        met_demand_per_sector : dictionary with sector names (string) as keys and PCRaster maps with
                                sectoral met demands as values (units: m3/d)
        message_str           : output message
        '''
        
        # [ set up input data ] .............................................................................................
        # define starting conditions
        remaining_availability      = deepcopy(availability)
        unmet_demand_per_sector     = deepcopy(demand_per_sector)
        
        met_demand_per_sector       = dict((sector_name, \
                                            pcr.ifthen(demand_per_sector[sector_name] >= 0.0, \
                                                       pcr.scalar(0.0))) \
                                           for sector_name in self.sector_names)
        
        withdrawal_per_sector       = dict((source_name, \
                                            dict((sector_name, \
                                                  pcr.spatial(pcr.scalar(0.0))) \
                                                 for sector_name in self.sector_names)) \
                                           for source_name in self.source_names)
        
        allocated_demand_per_sector = dict((source_name, \
                                            dict((sector_name, \
                                                  pcr.spatial(pcr.scalar(0.0))) \
                                                 for sector_name in self.sector_names)) \
                                           for source_name in self.source_names)
        
        # define allocation zones
        # domestic sector typically has the largest zone
        zones = {'surfacewater': zones_per_sector['surfacewater']['domestic'], \
                 'groundwater' : zones_per_sector['groundwater']['domestic']}
        
        # get long-term quality state for current date
        constituent_longterm_states = self.water_quality.get_longterm_quality_for_date( \
                                             source_names = self.source_names, \
                                             date         = date)
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        if verbose:
            dt = f'{str(date.year)[2:]}-{str(date.month).zfill(2)}'
            for source_name in self.source_names:
                for constituent_name in self.water_quality.constituent_names:
                    pcr.report(constituent_longterm_states[source_name][constituent_name], f'{path}/{dt}_longterm_{source_name}_{constituent_name}.map')
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # define suitability masks of water quality per sector and source
        self.suitability_per_sector = {}
        for source_name in self.source_names:
            self.suitability_per_sector[source_name] = \
                 self.water_quality.get_suitability_per_sector( \
                             constituent_state = constituent_longterm_states[source_name], \
                             sector_names      = self.sector_names)
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        if verbose:
            for source_name in self.source_names:
                for sector_name in self.sector_names:
                    pcr.report(self.suitability_per_sector[source_name][sector_name], f'{path}/{dt}_longterm_suitability_{source_name}_{sector_name}.map')
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # define initial exit conditions
        iter_allocation = 1
        exit_condition  = False
        max_iter_allocation = 25
        totz_demand_old = {}
        totz_supply_old = {}
        
        # [ start water distribution ] ......................................................................................
        # iterate until either the demand is met or the supply is exhausted
        while not exit_condition:
            
            # [ DELETEME ] verbose <------------------------------------------------------------------------------------------------------------------
            if verbose:
                dt = str(date)[2:7]
                itera = str(iter_allocation).zfill(2)
                print('date: %s - iteration: %s \n total available : surfacewater = %15s - groundwater  = %15s' % \
                      (dt, \
                       itera, \
                       pcr.cellvalue(pcr.maptotal(remaining_availability['surfacewater']),1)[0], \
                       pcr.cellvalue(pcr.maptotal(remaining_availability['groundwater']),1)[0]))
            # ----------------------------------------------------------------------------------------------------------------------------------------
            
            # define initial values of zonal demand and supply
            for source_name in self.source_names:
                if iter_allocation == 1:
                    totz_demand_old[source_name] = get_zonal_total( \
                                                       local_values = pcr.max(0, \
                                                                              sum_list(list(unmet_demand_per_sector.values()))), \
                                                       zones        = zones[source_name])
                    
                    totz_supply_old[source_name] = get_zonal_total( \
                                                       local_values = remaining_availability[source_name], \
                                                       zones        = zones[source_name])
                else:
                    totz_demand_old[source_name] = deepcopy(total_zonal_demand[source_name])
                    totz_supply_old[source_name] = deepcopy(total_zonal_supply[source_name])
            
            # evaluate potential water withdrawal
            # define water availability weights per sector accounting for water quality
            
            # [ surfacewater ]
            source_name = 'surfacewater'
            unmet_demand_per_sector_surfacewater = deepcopy(unmet_demand_per_sector)
            weights_surfacewater_per_sector = \
                self.water_quality.get_weights_availability_per_sector( \
                         source_name               = source_name, \
                         sector_names              = self.sector_names, \
                         prioritization_per_sector = prioritization[source_name], \
                         suitability_per_sector    = self.suitability_per_sector[source_name], \
                         demand_per_sector         = unmet_demand_per_sector_surfacewater, \
                         availability              = remaining_availability[source_name], \
                         zones_per_sector          = zones_per_sector[source_name], \
                         )
            
            # [ groundwater ]
            source_name = 'groundwater'
            unmet_demand_per_sector_groundwater = deepcopy(unmet_demand_per_sector)
            
            # if thermoelectric sector or environmental flow requirements are evaluated
            for sector_name in self.sectors_local_surfacewater:
                if sector_name in self.sector_names:
                    unmet_demand_per_sector_groundwater[sector_name] = pcr.spatial(pcr.scalar(0))
            
            # get the weights
            weights_groundwater_per_sector = \
                self.water_quality.get_weights_availability_per_sector( \
                         source_name               = source_name, \
                         sector_names              = self.sector_names, \
                         prioritization_per_sector = prioritization[source_name], \
                         suitability_per_sector    = self.suitability_per_sector[source_name], \
                         demand_per_sector         = unmet_demand_per_sector_groundwater, \
                         availability              = remaining_availability[source_name], \
                         zones_per_sector          = zones_per_sector[source_name], \
                         )
            
            # allocate the sectoral gross demand for current period to the demand
            # and return the withdrawal, allocated demand and the met and unmet demands per sector
            for sector_name in self.sector_names:
                
                # calculate the assigned available water per sector
                remaining_availability_surfacewater_sector = \
                    remaining_availability['surfacewater'] * weights_surfacewater_per_sector[sector_name] * self.suitability_per_sector['surfacewater'][sector_name]
                
                remaining_availability_groundwater_sector = \
                    remaining_availability['groundwater']  * weights_groundwater_per_sector[sector_name]  * self.suitability_per_sector['groundwater'][sector_name]
                
                # [ DELETEME ] verbose <--------------------------------------------------------------------------------------------------------------
                if verbose:
                    pcr.report(remaining_availability_surfacewater_sector, f'{path}/{dt}_longterm_surfacewater_{sector_name}_iter{itera}.map')
                    pcr.report(remaining_availability_groundwater_sector,  f'{path}/{dt}_longterm_groundwater_{sector_name}_iter{itera}.map')
                # ------------------------------------------------------------------------------------------------------------------------------------
                
                # calculate water withdrawal and allocation
                tmp_withdrawal, tmp_allocated_demand, tmp_met_demand, tmp_unmet_demand, message_str = \
                    allocate_demand_to_availability_with_options( \
                        demand             = unmet_demand_per_sector[sector_name], \
                        availability       = {'surfacewater' : remaining_availability_surfacewater_sector, \
                                              'groundwater'  : remaining_availability_groundwater_sector}, \
                        zones              = {'surfacewater' : zones_per_sector['surfacewater'][sector_name], \
                                              'groundwater'  : zones_per_sector['groundwater'][sector_name]}, \
                        source_names       = self.source_names, \
                        use_local_first     = use_local_first, \
                        reallocate_surplus = reallocate_surplus)
                
                # [ DELETEME ] verbose <--------------------------------------------------------------------------------------------------------------
                if verbose:
                    print(' %s\n  - available   : surfacewater = %15s - groundwater  = %15s \n  - withdrawals : surfacewater = %15s - groundwater  = %15s \n  - met demand  : %s \n  - unmet demand: %s' % \
                          (sector_name, \
                           pcr.cellvalue(pcr.maptotal(remaining_availability_surfacewater_sector),1)[0], \
                           pcr.cellvalue(pcr.maptotal(remaining_availability_groundwater_sector),1)[0], \
                           pcr.cellvalue(pcr.maptotal(tmp_withdrawal['surfacewater']),1)[0], \
                           pcr.cellvalue(pcr.maptotal(tmp_withdrawal['groundwater']),1)[0], \
                           pcr.cellvalue(pcr.maptotal(tmp_met_demand),1)[0], \
                           pcr.cellvalue(pcr.maptotal(tmp_unmet_demand),1)[0]))
                # ------------------------------------------------------------------------------------------------------------------------------------
                
                # update water withdrawal and demand values
                met_demand_per_sector[sector_name]   = \
                            met_demand_per_sector[sector_name] + tmp_met_demand
                
                unmet_demand_per_sector[sector_name] = \
                    pcr.max(0,
                            unmet_demand_per_sector[sector_name] - tmp_met_demand)
                
                for source_name in self.source_names:
                    withdrawal_per_sector[source_name][sector_name] = \
                        withdrawal_per_sector[source_name][sector_name] + tmp_withdrawal[source_name]
                    
                    allocated_demand_per_sector[source_name][sector_name] = \
                        allocated_demand_per_sector[source_name][sector_name] + tmp_allocated_demand[source_name]
            
            # aggregate withdrawals and allocations by source
            withdrawal       = dict((source_name, \
                                     sum_list(list(withdrawal_per_sector[source_name].values()))) \
                                    for source_name in self.source_names)
            allocated_demand = dict((source_name, \
                                     sum_list(list(allocated_demand_per_sector[source_name].values()))) \
                                    for source_name in self.source_names)
            
            # update remaining water available
            for source_name in self.source_names:
                remaining_availability[source_name] = \
                    pcr.max(0, \
                            availability[source_name] - withdrawal[source_name])
            
            # update iteration number and exit condition
            total_zonal_demand = {}
            total_zonal_supply = {}
            for source_name in self.source_names:
                total_zonal_demand[source_name] = get_zonal_total( \
                                                      local_values = sum_list(list(unmet_demand_per_sector.values())), \
                                                      zones        = zones[source_name])
                total_zonal_supply[source_name] = get_zonal_total( \
                                                      local_values = remaining_availability[source_name], \
                                                      zones        = zones[source_name])
            
            update_mask     = ((total_zonal_demand['surfacewater'] < totz_demand_old['surfacewater']) & \
                               (total_zonal_supply['surfacewater'] < totz_supply_old['surfacewater'])) | \
                              ((total_zonal_demand['groundwater']  < totz_demand_old['groundwater']) & \
                               (total_zonal_supply['groundwater']  < totz_supply_old['groundwater']))
            
            iter_allocation = iter_allocation + 1
            exit_condition  = (pcr.cellvalue(pcr.mapmaximum(pcr.scalar(update_mask)), 1)[0] == 0) | \
                              (iter_allocation > max_iter_allocation)
            
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
            if verbose:
                print('  exit condition -> \n   - surfacewater : demands = %s - availability = %s -> overall = %s\n   - groundwater  : demands = %s - availability = %s -> overall = %s\n   - final overall : %s' % \
                      (pcr.cellvalue(pcr.mapmaximum(pcr.scalar((total_zonal_demand['surfacewater'] < totz_demand_old['surfacewater']))), 1)[0] == 0, \
                       pcr.cellvalue(pcr.mapmaximum(pcr.scalar((total_zonal_supply['surfacewater'] < totz_supply_old['surfacewater']))), 1)[0] == 0, \
                       pcr.cellvalue(pcr.mapmaximum(pcr.scalar((total_zonal_demand['surfacewater'] < totz_demand_old['surfacewater']) & (total_zonal_supply['surfacewater'] < totz_supply_old['surfacewater']))), 1)[0] == 0, \
                       pcr.cellvalue(pcr.mapmaximum(pcr.scalar((total_zonal_demand['groundwater']  < totz_demand_old['groundwater']))), 1)[0] == 0, \
                       pcr.cellvalue(pcr.mapmaximum(pcr.scalar((total_zonal_supply['groundwater']  < totz_supply_old['groundwater']))), 1)[0] == 0, \
                       pcr.cellvalue(pcr.mapmaximum(pcr.scalar((total_zonal_demand['groundwater']  < totz_demand_old['groundwater']) & (total_zonal_supply['groundwater']  < totz_supply_old['groundwater']))), 1)[0] == 0, \
                       exit_condition))
        
        if verbose:
            print(' outcomes \n - withdrawals : surfacewater = %15s - groundwater  = %15s \n - met demand  : %s \n - unmet demand: %s' % \
                  (pcr.cellvalue(pcr.maptotal(withdrawal['surfacewater']),1)[0], \
                   pcr.cellvalue(pcr.maptotal(withdrawal['groundwater']),1)[0], \
                   pcr.cellvalue(pcr.maptotal(sum_list(list(met_demand_per_sector.values()))),1)[0], \
                   pcr.cellvalue(pcr.maptotal(sum_list(list(unmet_demand_per_sector.values()))),1)[0]))
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # return potential withdrawals and met demand per sector
        return withdrawal_per_sector, met_demand_per_sector, message_str
    
    
    
    def allocate_unmet_demand_to_nonrenewable_sources(self, \
                                                       unmet_demand_per_sector, \
                                                       zones_per_sector, \
                                                       withdrawal_capacity, \
                                                       withdrawal_points, \
                                                       ):
        '''
        allocate_unmet_demand_to_nonrenewable_sources:
                                       function to allocate the unmet demand (or claimed withdrawal)
                                       per sector to the groundwater sources based on its potential
                                       renewable withdrawals
        
        input:
        =====
        unmet_demand_per_sector      : dictionary with sector names (string) as keys and PCRaster maps
                                       with the amount of withdrawal that cannot be met and are assigned
                                       to the non-renewable groundwater (units: m3/day)
        zones_per_sector             : dictionary with sector names (string) as keys and PCRaster maps
                                       with allocation zones per sector (nominal)
        withdrawal_capacity          : dictionary with source names (string) as keys and PCRaster maps
                                       with water withdrawal capacity as values (units: m3/day)
        withdrawal_points            : dictionary with source names (string) as keys and PCRaster
                                       maps with water withdrawal points as values (nominal)
        '''
        
        message_str = 'Non-renewable withdrawals are assigned to the following sources'
        
        # set variables
        source_name         = 'groundwater'
        availability        = deepcopy(self.potential_renewable_withdrawal[source_name])
        suitability         = deepcopy(self.suitability_per_sector[source_name])
        zones               = deepcopy(zones_per_sector[source_name])
        withdrawal_capacity = deepcopy(withdrawal_capacity[source_name])
        
        # filter sectors that can withdraw water also from groundwater
        # and unmet demands from these sectors
        sector_names = deepcopy(self.sector_names)
        for sector_name in self.sectors_local_surfacewater:
            if sector_name in self.sector_names:
                sector_names.remove(sector_name)
                unmet_demand_per_sector.pop(sector_name)
        
        # define remaining withdrawal capacity (if defined) per sector
        # note:
        # potential non-renewable withdrawal is also accounted for,
        # being this variable zero during the beginning of the allocation process
        if not isinstance(withdrawal_capacity, NoneType):
            withdrawal_capacity_remaining = \
                pcr.max(0, \
                        withdrawal_capacity - (self.potential_renewable_withdrawal[source_name] + \
                                               self.potential_nonrenewable_withdrawal[source_name]))
            
            withdrawal_capacity_remaining_per_sector = {}
            for sector_name in sector_names:
                sector_rate     = pcr_return_val_div_zero( \
                                          get_zonal_total(unmet_demand_per_sector[sector_name], \
                                                          zones[sector_name]),
                                          get_zonal_total(sum_list(list(unmet_demand_per_sector.values())), \
                                                          zones[sector_name]),
                                          very_small_number)
                
                withdrawal_capacity_remaining_per_sector[sector_name] = \
                                      withdrawal_capacity_remaining * sector_rate
        
        # get the availability per sector considering the water quality,
        # note that this favours pumping capacity in case water is unavailable
        # so highlights non-renewable withdrawals
        for sector_name in sector_names:
            
            # define the withdrawal capacity
            # if it is not defined, it takes the minimum availability value if
            # it is not zero, else it makes it 1 to ensure water can always be
            # withdrawn
            if not isinstance(withdrawal_capacity, NoneType):
                capacity = withdrawal_capacity_remaining
                availability_per_sector = capacity * suitability[sector_name]
            else:
                capacity = pcr.areaminimum(availability, zones[sector_name])
                capacity = pcr.ifthenelse(capacity > 0, capacity, pcr.scalar(1))
                
                availability_per_sector = pcr.ifthenelse(availability > 0, \
                                                         availability * suitability[sector_name], \
                                                         capacity     * suitability[sector_name])
            
            # define the availability considering the sectoral water quality
            # requirements only from the points where water can be withdrawn
            availability_per_sector = \
                pcr.ifthenelse( \
                               pcr.scalar(withdrawal_points[source_name]) > 0, \
                               availability_per_sector, \
                               0)
            
            # get the potential non-renewable withdrawal per sector
            # (units: m3/day)
            allocation_rate = pcr_return_val_div_zero( \
                                              availability_per_sector, \
                                              get_zonal_total(availability_per_sector, \
                                                              zones[sector_name]), \
                                              very_small_number)
            
            allocated_demand = allocation_rate * \
                               get_zonal_total(unmet_demand_per_sector[sector_name], \
                                               zones[sector_name])
            
            # evaluate if pumping capacity limits the allocated demand
            if not isinstance(withdrawal_capacity, NoneType):
                allocated_demand = pcr.min(allocated_demand, \
                                           withdrawal_capacity_remaining_per_sector[sector_name])
            
            # set variable (units: m3/day)
            self.potential_nonrenewable_withdrawal_per_sector[source_name][sector_name] = \
                              self.potential_nonrenewable_withdrawal_per_sector[source_name][sector_name] + \
                              allocated_demand
        
        # aggregate potential non-renewable withdrawal from all sectors
        # (units: m3/day)
        self.potential_nonrenewable_withdrawal[source_name] = \
             sum_list(list(self.potential_nonrenewable_withdrawal_per_sector[source_name].values()))
        
        # log message
        message_str = str.join(': ', \
                               (message_str, source_name))
        logger.debug(message_str)
        
        # returns None
        return None
    
    
    
    def get_potential_withdrawal(self, \
                                 source_name):
        '''
        get_potential_withdrawal:
                                  function to get the potential withdrawals per sector
                                  as the sum of renewable and non-renewable water withdrawals
                                  limited by the withdrawal capacity
        
        input:
        =====
        source_name             : string, source name under evaluation
        
        output:
        ======
        total_potential_withdrawal_per_sector:
                                  dictionary with sector names (string) as keys and
                                  PCRaster maps with sum of renewable and non-renewable
                                  potential withdrawals per sector, limited by withdrawal
                                  capacity (units: m3/day)
        '''
        
        message_str = 'updating potential withdrawals per sector'
        
        # verify potential withdrawals are affected by withdrawal capacity
        if not isinstance(self.surfacewater_withdrawal_capacity, NoneType):
            if pcr.cellvalue(pcr.mapminimum(self.surfacewater_withdrawal_capacity - \
                                            (self.potential_renewable_withdrawal['surfacewater'] + \
                                             self.potential_nonrenewable_withdrawal['surfacewater'])), 1)[0] < -1:
                logger.info('Sum of potential surface water renewable and non-renewable withdrawals are larger than surface water withdrawal capacity')
                sys.exit()
        
        if not isinstance(self.groundwater_withdrawal_capacity, NoneType):
            if pcr.cellvalue(pcr.mapminimum(self.groundwater_withdrawal_capacity - \
                                            (self.potential_renewable_withdrawal['groundwater'] + \
                                             self.potential_nonrenewable_withdrawal['groundwater'])), 1)[0] < -1:
                logger.info('Sum of potential groundwater renewable and non-renewable withdrawals are larger than groundwater withdrawal capacity')
                sys.exit()
        
        # get total long-term potential renewable and non-renewable withdrawals per sector
        # (units: m3/day)
        total_potential_withdrawal_per_sector = \
                       dict((sector_name, \
                             self.potential_renewable_withdrawal_per_sector[source_name][sector_name] + \
                             self.potential_nonrenewable_withdrawal_per_sector[source_name][sector_name])
                            for sector_name in self.sector_names)
        
        # log message string
        logger.debug(message_str)
        
        # return the total potential withdrawal (units: m3/day)
        return total_potential_withdrawal_per_sector
    
    
    
    def update_surfacewater_potential_withdrawals(self, \
                                                   total_runoff, \
                                                   surfacewater_storage, \
                                                   fraction_water, \
                                                   longterm_potential_withdrawal_per_sector, \
                                                   cellarea, \
                                                   date = None):
        '''
        update_surfacewater_potential_withdrawals:
                                       function that calculates the actual water withdrawals from the
                                       channel based on the potential sectoral demands (potential_withdrawal)
                                       and the oustanding potential evapotranspiration (channel_runoff < 0)
        
        input:
        =====
        total_runoff                  : PCRaster map with available runoff as the sum of the different
                                       surface water components (units: m/day)
        surfacewater_storage         : PCRaster map with surface water storage at the start of the time-step
                                       (units: m per day)
        fraction_water               : PCRaster map with fraction of the area of a pixel covered by water
                                       (unitless)
        longterm_potential_withdrawal_per_sector :
                                       dictionary with sector names (string) as keys and PCRaster maps
                                       with sum of renewable and non-renewable potential withdrawals per
                                       sector obtained considering long-term water quality as values
                                       (units: m3/day)
        cellarea                     : PCRaster map with cell area (units: m2)
        date                         : string, date under evaluation
        
        output:
        ======
        actual_withdrawal_per_sector : dictionary with sector names (string) as keys and PCRaster maps 
                                       with actual water withdrawal from sectoral demands based on water
                                       availability (units: m/period)
        channel_runoff_per_sector     : dictionary with sector names (string) as keys and PCRaster maps 
                                       with channel runoff with actual evapotranspiration from channel based
                                       on water availability (units: m/period)
        '''
        
        source_name = 'surfacewater'
        
        # get suitability per sector considering short-term surface water quality
        # (units: -)
        suitability_per_sector = self.water_quality.get_suitability_per_sector( \
                 constituent_state = self.water_quality.constituent_shortterm_quality[source_name], \
                 sector_names      = self.sector_names)
        
        # get the short-term potential surface water availability
        # (units: m3/day)
        potential_surfacewater_availability = \
                 (surfacewater_storage * fraction_water + total_runoff ) * cellarea
        
        # get the short-term potential surface water withdrawals
        # by updating the long-term potential surface water withdrawals 
        # considering the short-term water quality suitability
        # (units: m3/day)
        shortterm_potential_withdrawal_per_sector = \
            dict((sector_name, \
                  longterm_potential_withdrawal_per_sector[sector_name] * \
                   suitability_per_sector[sector_name]) \
                 for sector_name in self.sector_names)
        
        # aggregate potential withdrawals from all sectors
        # (units: m3/day)
        shortterm_potential_withdrawal = \
                 sum_list(list(shortterm_potential_withdrawal_per_sector.values()))
        
        # get the current potential surface water withdrawals based on the
        # suitability per sector considering the short-term quality
        # (units: m3/day)
        potential_withdrawal_per_sector = \
            dict((sector_name, \
                  pcr.min(shortterm_potential_withdrawal_per_sector[sector_name], \
                          potential_surfacewater_availability * \
                           suitability_per_sector[sector_name] * \
                           pcr_return_val_div_zero(shortterm_potential_withdrawal_per_sector[sector_name], \
                                                   shortterm_potential_withdrawal, \
                                                   very_small_number))) \
                 for sector_name in self.sector_names)
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        if verbose:
            dt = f'{str(date.year)[2:]}-{str(date.month).zfill(2)}'
            pcr.report(potential_surfacewater_availability, f'{path}/{dt}_shortterm_potential_availability_surfacewater.map')
            for sector_name in self.sector_names:
                pcr.report(suitability_per_sector[sector_name], f'{path}/{dt}_shortterm_suitability_surfacewater_{sector_name}.map')
                pcr.report(potential_withdrawal_per_sector[sector_name], f'{path}/{dt}_shortterm_potential_withdrawal_surfacewater_{sector_name}.map')
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # return potential surface water withdrawals per sector
        # for current time-step
        # (units: m3/day)
        return potential_withdrawal_per_sector
    
    
    
    def update_groundwater_actual_withdrawals(self, \
                                              groundwater_storage, \
                                              total_recharge, \
                                              total_base_flow, \
                                              longterm_potential_withdrawal_per_sector, \
                                              cellarea, \
                                              time_step_length, \
                                              date=None):
        '''
        update_groundwater_actual_withdrawals:
                                             function to update the potential groundwater withdrawals
                                             considering the short-term water quality and distributing
                                             among renewable and non-renewable sources based on the
                                             storage
        
        input:
        =====
        groundwater_storage                : PCRaster map with groundwater storage of the previous time-
                                             step (units: m)
        total_recharge                     : PCRaster map with groundwater recharge accumulated at the
                                             end of the period (units: m/period)
        total_base_flow                     : PCRaster map with groundwater base flow accumulated at the
                                             end of the period (units: m/period)
        longterm_potential_withdrawal_per_sector : 
                                             dictionary with sector names (string) as keys and PCRaster maps
                                             with sum of renewable and non-renewable potential withdrawals per
                                             sector obtained considering long-term water quality as values
                                             (units: m/period)
        cellarea                           : PCRaster map with area of cells (units: m2)
        time_step_length                   : integer, number of days in period (e.g., 30 days/month)
        
        output:
        ======
        storage                            : PCRaster map with groundwater renewable storage before 
                                             water withdrawals (units m)
        renewable_withdrawal_per_sector    : PCRaster map with renewable withdrawals per sector
                                             (units: m3/day)
        nonrenewable_withdrawal_per_sector : PCRaster map with non-renewable withdrawals per sector
                                             (units: m3/day)
        '''
        source_name = 'groundwater'
        
        # [ groundwater availability ]
        # update the storage at the beginning of the time-step
        # with the total recharge and total base flow at the end of the time-step
        # (units: m at the end of the period)
        storage = groundwater_storage + total_recharge - total_base_flow
        
        # get the short-term groundwater availability
        # (units: m3 at the end of the period)
        groundwater_availability = pcr.max(0, storage * cellarea)
        
        # [ short-term potential withdrawal ]
        # get suitability per sector considering short-term groundwater quality
        suitability_per_sector = \
            self.water_quality.get_suitability_per_sector( \
                     constituent_state = self.water_quality.constituent_shortterm_quality[source_name], \
                     sector_names      = self.sector_names)
        
        # update long-term potential groundwater withdrawals considering 
        # short-term water quality suitability
        # (units: m3/period)
        potential_withdrawal_per_sector = \
            dict((sector_name, \
                  longterm_potential_withdrawal_per_sector[sector_name] \
                  * suitability_per_sector[sector_name] \
                  * time_step_length) \
                 for sector_name in self.sector_names)
        
        # aggregate total potential withdrawals from all sectors
        # (units: m3/period)
        potential_withdrawal = \
                 sum_list(list(potential_withdrawal_per_sector.values()))
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        if verbose:
            dt = f'{str(date.year)[2:]}-{str(date.month).zfill(2)}'
            pcr.report(groundwater_availability / time_step_length, f'{path}/{dt}_shortterm_actual_renewable_availability_groundwater.map')
            pcr.report(potential_withdrawal / time_step_length, f'{path}/{dt}_shortterm_potential_withdrawal_groundwater.map')
            for sector_name in self.sector_names:
                pcr.report(suitability_per_sector[sector_name], f'{path}/{dt}_shortterm_suitability_groundwater_{sector_name}.map')
                pcr.report(potential_withdrawal_per_sector[sector_name] / time_step_length, f'{path}/{dt}_shortterm_potential_withdrawal_groundwater_{sector_name}.map')
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # [ actual withdrawals ]
        # calculate the renewable and non-renewable withdrawals
        # (units: m3/period)
        renewable_withdrawal = \
                  pcr.min(groundwater_availability, \
                          potential_withdrawal)
        
        nonrenewable_withdrawal = \
                  pcr.max(0, \
                          potential_withdrawal - renewable_withdrawal)
        
        # convert renewable and non-renewable withdrawals
        # (units: m3/day)
        renewable_withdrawal    = renewable_withdrawal / time_step_length
        nonrenewable_withdrawal = nonrenewable_withdrawal / time_step_length
        
        # re-distribute renewable withdrawals by sector
        # (units: m3/day)
        renewable_withdrawal_per_sector = \
             dict((sector_name, \
                   renewable_withdrawal * \
                    pcr_return_val_div_zero(potential_withdrawal_per_sector[sector_name], \
                                            potential_withdrawal, \
                                            very_small_number)) \
                  for sector_name in self.sector_names)
        
        # re-distribute non-renewable withdrawals by sector (units: m3/day)
        nonrenewable_withdrawal_per_sector = \
             dict((sector_name, \
                   nonrenewable_withdrawal * \
                    pcr_return_val_div_zero(potential_withdrawal_per_sector[sector_name], \
                                            potential_withdrawal, \
                                            very_small_number)) \
                  for sector_name in self.sector_names)
        
        # return actual withdrawals from renewable and non-renewable sources,
        # per sector (units: m3/day) and renewable groundwater storage (units: m)
        return renewable_withdrawal_per_sector, nonrenewable_withdrawal_per_sector, \
               renewable_withdrawal, nonrenewable_withdrawal, storage
                        
    
    
    
    def update_withdrawals(self, \
                            source_name, \
                            renewable_withdrawal_per_sector, \
                            nonrenewable_withdrawal_per_sector, \
                            source_names_to_be_processed, \
                            date=None):
        '''
        update_withdrawals           : function that wraps around two individual actions
                                       that are needed to update the potential and actual
                                       withdrawals:
                                       1) set the actual withdrawals on the basis of the
                                          water that can be withdrawn from the source provided
                                          in a renewable or non-renewable fashion;
                                       2) pass any unmet demand to the nonrenewable potential
                                          withdrawals to the remaining allowable sources.
        
        input:
        =====
        source_name                  : source name currently processed; the renewable
                                       and non-renewable withdrawals are the amounts
                                       that currently could be actually withdrawn
        renewable_withdrawal_per_sector,
        nonrenewable_withdrawal_per_sector : 
                                       total amount of water withdrawn from the
                                       renewable and non-renewable storage for the
                                       current source (units: m3/period)
        source_names_to_be_processed : a list of the sector names that are elligi-
                                       ble to accommodate any unmet demand as pot-
                                       ential non-renewable withdrawals.
        '''
        
        # [ actual withdrawals ]
        self.set_actual_withdrawals(source_name, \
                                    renewable_withdrawal_per_sector, \
                                    nonrenewable_withdrawal_per_sector)
        
        # [ unmet demand ]
        # set the unmet demand per sector as the difference of what was potentially 
        # withdrawn and what is actually withdrawn
        unmet_withdrawal_per_sector = \
            dict((sector_name,
                  pcr.max(0, \
                          (self.potential_renewable_withdrawal_per_sector[source_name][sector_name] + \
                           self.potential_nonrenewable_withdrawal_per_sector[source_name][sector_name]) \
                            - \
                          (self.actual_renewable_withdrawal_per_sector[source_name][sector_name] + \
                           self.actual_nonrenewable_withdrawal_per_sector[source_name][sector_name])) )\
                 for sector_name in self.sector_names)
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        if verbose:
            dt = f'{str(date.year)[2:]}-{str(date.month).zfill(2)}'
            for sector_name in self.sector_names:
                pcr.report(unmet_withdrawal_per_sector[sector_name], f'{path}/{dt}_unmet_{source_name}_withdrawal_{sector_name}_after_{source_name}_[2reallocate].map')
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # set a list of selectable source names to process
        source_names = []
        for source_name in source_names_to_be_processed:
            source_names.append(source_name)
        
        # [ update non-renewable withdrawals ]
        # allocate unmet withdrawals to non-renewable groundwater source
        if 'groundwater' in source_names:
            message_str = 'Non-renewable withdrawals are assigned to the following sources: groundwater' \
            
            self.allocate_unmet_demand_to_nonrenewable_sources( \
                     unmet_demand_per_sector = unmet_withdrawal_per_sector,\
                     zones_per_sector        = {'surfacewater' : self.surfacewater_allocation_zones, \
                                                'groundwater'  : self.groundwater_allocation_zones}, \
                     withdrawal_capacity     = {'surfacewater' : self.surfacewater_withdrawal_capacity, \
                                                'groundwater'  : self.groundwater_withdrawal_capacity}, \
                     withdrawal_points       = {'surfacewater' : self.surfacewater_withdrawal_points, \
                                                'groundwater'  : self.groundwater_withdrawal_points})
            # log message
            logger.debug(message_str)
        
        # returns None
        return None
    
    
    
    def set_actual_withdrawals(self, \
                                source_name, \
                                renewable_withdrawal_per_sector, \
                                nonrenewable_withdrawal_per_sector):
        
        '''
        set_actual_withdrawals             : function that sets the actual withdrawals for
                                             the source given the name and the amounts of
                                             renewable and non-renewable withdrawals in a cell
        
        input:
        =====
        source_name                        : name of the source currently being processed;
                                             for this source the renewable and non-renewable
                                             withdrawals are the amounts thar currently could
                                             be actually withdrawn
        renewable_withdrawal_per_sector    : total amount of water withdrawn from the
                                             renewable storage for the current source;
        nonrenewable_withdrawal_per_sector : idem, but then withdrawn from the non-
                                             renewable storage for the current source.
        '''

        # function to set the actual withdrawals
        self.actual_renewable_withdrawal[source_name]    = \
                 sum_list(list(renewable_withdrawal_per_sector.values()))
        
        self.actual_nonrenewable_withdrawal[source_name] = \
                 sum_list(list(nonrenewable_withdrawal_per_sector.values()))
        
        self.actual_renewable_withdrawal_per_sector[source_name] = \
                 dict((sector_name, \
                       renewable_withdrawal_per_sector[sector_name]) \
                      for sector_name in self.sector_names)
        
        self.actual_nonrenewable_withdrawal_per_sector[source_name] = \
                 dict((sector_name, \
                       nonrenewable_withdrawal_per_sector[sector_name]) \
                      for sector_name in self.sector_names)
        
        # returns None
        return None
    
    
    
    def allocate_withdrawal_to_demand_for_date(self, date):
        
        '''
        allocate_demand_to_withdrawals:
               function that internally allocates the gross demand to the sources
               on a cell-by-cell basis given the actual withdrawals and updates
               the consumption and return flows  
        
        input:
        =====
        date : string, date of the update
        '''
        # set the message string to log the information
        message_str = 'Actual water withdrawals allocated to the demand for %s.' % date
        
        # iterate over the renewable and non-renewable withdrawals and per source
        # to allocate the water demand; allocated withdrawal and demand are set
        # as internal variables, unused withdrawal and met demand are not used
        # directly; the remaining, unused withdrawals are kept in a structure
        # similar to the potential and actual demand for checks, however, and
        # added to the return flows to avoid balance errors; the sub_message_str
        # is the text string with the information on the allocation for inspec-
        # tion and logging (units: m3/day)
        self.allocated_withdrawal_per_sector, unused_withdrawal, \
        self.allocated_demand_per_sector,  met_demands, sub_message_str = \
               allocate_demand_to_withdrawals( \
                    withdrawal_names                   = self.withdrawal_names, \
                    source_names                       = self.source_names, \
                    sector_names                       = self.sector_names, \
                    demand_per_sector                  = self.gross_demand_remaining, \
                    renewable_withdrawal_per_sector    = self.actual_renewable_withdrawal_per_sector, \
                    nonrenewable_withdrawal_per_sector = self.actual_nonrenewable_withdrawal_per_sector, \
                    zones_per_sector                   = {'surfacewater' : self.surfacewater_allocation_zones, \
                                                          'groundwater'  : self.groundwater_allocation_zones}, \
                    use_local_first                     = self.use_local_first, \
                    )
        
        # update the message_str
        message_str = str.join('\n', \
                               (message_str, sub_message_str))
        
        # NOTE: this is a bit silly but just to keep things tractable:
        # add the unused withdrawals per withdrawal type and source
        # idem for the actual ones (units: m3/day)
        for withdrawal_name in self.withdrawal_names:
                
            var_str = get_key(['unused', withdrawal_name, 'withdrawal'])
            
            setattr(self, var_str, dict((source_name, \
                                        unused_withdrawal[withdrawal_name][source_name]) \
                                   for source_name in self.source_names))
        
        # iterate over the information on the return flows and the consumption
        # and the withdrawal and alllocations
        # total withdrawal is dependent on the actual withdrawals which also
        # includes unused withdrawals (units: m3/day)
        self.total_withdrawal = sum_list(list(self.actual_renewable_withdrawal.values())) + \
                                sum_list(list(self.actual_nonrenewable_withdrawal.values()))
        
        # initialize the total allocated water, consumption and the return flow
        # that are updated by iterating over the sectors
        self.total_allocation  = pcr.spatial(pcr.scalar(0))
        self.total_consumption = pcr.spatial(pcr.scalar(0))
        self.total_return_flow = pcr.spatial(pcr.scalar(0))
        
        for sector_name in self.sector_names:
            
            # get the return flow ratio (units: -)
            return_flow_ratio = self.get_return_flow_ratio( \
                                         gross_demand = self.gross_demand[sector_name], \
                                         net_demand   = self.net_demand[sector_name])
            
            # iterate over the keys in the allocated demand (units: m3/day)
            for key in self.allocated_demand_per_sector.keys():
                
                # update the consumption and return flow per sector (units: m3/day)
                self.return_flow_demand_per_sector[key][sector_name] = \
                                 return_flow_ratio  * self.allocated_demand_per_sector[key][sector_name]
                
                self.consumed_demand_per_sector[key][sector_name] = \
                                 pcr.max(0, \
                                         self.allocated_demand_per_sector[key][sector_name] - \
                                         self.return_flow_demand_per_sector[key][sector_name])
                
                # update the total allocation (units: m3/day)
                self.total_allocation = self.total_allocation + \
                                        self.allocated_demand_per_sector[key][sector_name]
                    
                # update the total consumption and return flow (units: m3/day)
                self.total_return_flow = self.total_return_flow + \
                                         self.return_flow_demand_per_sector[key][sector_name]
                self.total_consumption = self.total_consumption + \
                                         self.consumed_demand_per_sector[key][sector_name]
            
            # include water use from desalinated water source if used
            if self.desalinated_water_use_flag:
                self.return_flow_demand_per_sector_desalwater[sector_name] = \
                                 return_flow_ratio  * self.allocated_demand_per_sector_desalwater[sector_name]
                self.consumed_demand_per_sector_desalwater[sector_name] = pcr.max(0, \
                                 self.allocated_demand_per_sector_desalwater[sector_name] - \
                                 self.return_flow_demand_per_sector_desalwater[sector_name])
                
                # update the total allocation
                self.total_allocation = self.total_allocation + \
                                        self.allocated_demand_per_sector_desalwater[sector_name]
                    
                # update the total consumption and return flow
                self.total_return_flow = self.total_return_flow + \
                                         self.return_flow_demand_per_sector_desalwater[sector_name]
                self.total_consumption = self.total_consumption + \
                                         self.consumed_demand_per_sector_desalwater[sector_name]
        
        # total return flow also contains the unused withdrawals (units: m3/day)
        self.total_return_flow = self.total_return_flow + \
            sum_list(list(self.unused_renewable_withdrawal.values())) + \
            sum_list(list(self.unused_nonrenewable_withdrawal.values()))
        
        # log the message
        logger.info(message_str)
        # add the final message
        logger.info('return flows and consumption added on the basis of the allocated demand')
        
        # returns None
        return None
    
    
    
    def get_return_flow_ratio(self, gross_demand, net_demand):
        '''
        get_return_flow_ratio: 
                           function which returns the return flow ratio as the
                           ratio of the net and gross demand per sector.
                           Assumes that all input is compatible with spatial scalar PCRaster fields.
        
        input:
        =====
        gross_demand     : gross water demand per sector [volume or waterslice]
        net_demand       : net water demand per sector   [volume or waterslice]
        
        output:
        ======
        return_flow_ratio : return flow ratio [-]
        '''
        
        # get the return flow ratio; small values of gross water demand result
        # in a return flow ratio of zero
        return_flow_ratio = pcr.max(0.00, 1.00 - \
                                    pcr_return_val_div_zero(net_demand, \
                                                            gross_demand, \
                                                            very_small_number))
        # return the return flow ratio
        return return_flow_ratio
    
    
    
    def update_longterm_availability(self, \
                                      groundwater_storage, \
                                      surfacewater_discharge, \
                                      surfacewater_runoff, \
                                      date):
        '''
        update_longterm_availability: function that updates the availability per
        zone as a function of the date.
        
        input:
        ======
        groundwater_storage    : PCRaster maps with groundwater potential storage
                                 the last fay of the month (units: m)
        surfacewater_discharge : PCRaster maps with surface discharge (units: m3/s)
        surfacewater_runoff     : PCRaster maps with surface total runoff (units: m/day)
        date                   : date of the update
        '''
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        dt = f'{str(date.year)[2:]}-{str(date.month).zfill(2)}'
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # [ groundwater storage ] ..................................................................
        # get the time step to update the groundwater storage
        date_index, matched_date, message_str = match_date_by_julian_number(date, \
                                                      self.groundwater_longterm_storage_dates)
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        if verbose:
            pcr.report(self.groundwater_longterm_storage[matched_date], f'{path}/{dt}_longterm_groundwater_storage_[initial].map')
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # remove the date from the dictionary and update it with the present value
        # set the value using the weight, if the long-term availability is not
        # defined, cover with the present value
        # (units: m per day)
        groundwater_longterm_storage = self.groundwater_longterm_storage.pop(matched_date)
        groundwater_longterm_storage = \
            pcr.cover( self.groundwater_update_weight      * groundwater_storage  + \
                      (1 - self.groundwater_update_weight) * groundwater_longterm_storage, \
                      groundwater_storage)
        # reset the date
        self.groundwater_longterm_storage_dates[date_index] = date
        
        # add the value to the dictionary
        self.groundwater_longterm_storage[date] = groundwater_longterm_storage
        
        # echo to screen
        message_str = str.join(' ', \
                                ('groundwater long-term total base flow updated for', \
                                message_str))
        logger.debug(message_str)
        
        # [ surface water discharge ] ..............................................................
        # get the time step to update the surface water availability (monthly)
        date_index, matched_date, message_str = match_date_by_julian_number(date, \
                                                      self.surfacewater_longterm_discharge_dates)
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        if verbose:
            pcr.report(self.surfacewater_longterm_discharge[matched_date], f'{path}/{dt}_longterm_surfacewater_discharge_[initial].map')
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # remove the date from the dictionary and update it with the present value
        # set the value using the weight, if the long-term availability is not
        # defined, cover with the present value
        # (units: m3/s)
        surfacewater_longterm_discharge = self.surfacewater_longterm_discharge.pop(matched_date)
        surfacewater_longterm_discharge = \
            pcr.cover( self.surfacewater_update_weight      * surfacewater_discharge + \
                      (1 - self.surfacewater_update_weight) * surfacewater_longterm_discharge, \
                      surfacewater_discharge)
        
        # reset the date
        self.surfacewater_longterm_discharge_dates[date_index] = date
        
        # add the value to the dictionary
        self.surfacewater_longterm_discharge[date] = surfacewater_longterm_discharge
        
        # echo to screen
        message_str = str.join(' ', \
                                ('surface water long-term discharge updated for', \
                                message_str))
        logger.debug(message_str)
        
        # [ surface water runoff ] .................................................................
        # get the time step to update the surface water availability (monthly)
        date_index, matched_date, message_str = match_date_by_julian_number(date, \
                                                      self.surfacewater_longterm_runoff_dates)
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        if verbose:
            pcr.report(self.surfacewater_longterm_runoff[matched_date], f'{path}/{dt}_longterm_surfacewater_runoff_[initial].map')
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # remove the date from the dictionary and update it with the present value
        # set the value using the weight, if the long-term availability is not
        # defined, cover with the present value
        # (units: m/day)
        surfacewater_longterm_runoff = self.surfacewater_longterm_runoff.pop(matched_date)
        surfacewater_longterm_runoff = \
            pcr.cover( self.surfacewater_update_weight      * surfacewater_runoff + \
                      (1 - self.surfacewater_update_weight) * surfacewater_longterm_runoff, \
                      surfacewater_runoff)
        
        # reset the date
        self.surfacewater_longterm_runoff_dates[date_index] = date
        
        # add the value to the dictionary
        self.surfacewater_longterm_runoff[date] = surfacewater_longterm_runoff
        
        # [ DELETEME ] verbose <----------------------------------------------------------------------------------------------------------------------
        if verbose:
            dt = f'{str(date.year)[2:]}-{str(date.month).zfill(2)}'
            pcr.report(groundwater_storage, f'{path}/{dt}_longterm_groundwater_storage_[final].map')
            pcr.report(surfacewater_discharge, f'{path}/{dt}_longterm_surfacewater_discharge_[final].map')
            pcr.report(surfacewater_runoff, f'{path}/{dt}_longterm_surfacewater_runoff_[final].map')
            pcr.report(groundwater_longterm_storage, f'{path}/{dt}_longterm_groundwater_storage_[updated].map')
            pcr.report(surfacewater_longterm_discharge, f'{path}/{dt}_longterm_surfacewater_discharge_[updated].map')
            pcr.report(surfacewater_longterm_runoff, f'{path}/{dt}_longterm_surfacewater_runoff_[updated].map')
        # --------------------------------------------------------------------------------------------------------------------------------------------
        
        # echo to screen
        message_str = str.join(' ', \
                                ('surface water long-term total runoff updated for', \
                                message_str))
        logger.debug(message_str)
        
        # returns None
        return None
    
    
    
    def update_longterm_potential_withdrawals(self, \
                                               date):
        '''
        update_longterm_potential_withdrawals: 
               function that updates the groundwater potential withdrawals as a function
               of the date
        
        input:
        ======
        date : date of the update.
        '''
        
        # [ groundwater ]
        #if self.pumping_capacity_flag['groundwater']:
        # get groundwater potential withdrawals for date
        # (units: m3/day)
        groundwater_potential_withdrawal = self.groundwater_potential_estimated_withdrawal
        
        # get the time step to update the groundwater potential withdrawal (monthly)
        # variable matches with groundwater_longterm_avail_dates
        date_index, matched_date, message_str = \
                    match_date_by_julian_number(date, \
                                                self.groundwater_longterm_pot_withdrawal_dates)
        
        # remove the date from the dictionary and update it with the present value
        # set the value using the weight, if the long-term availability is not
        # defined, cover with the present value
        groundwater_longterm_pot_withdrawal = self.groundwater_longterm_potential_withdrawal.pop(matched_date)
        groundwater_longterm_pot_withdrawal = \
              pcr.cover( self.groundwater_update_weight      * groundwater_potential_withdrawal + \
                        (1 - self.groundwater_update_weight) * groundwater_longterm_pot_withdrawal, \
                        groundwater_potential_withdrawal)
        
        # reset the date
        self.groundwater_longterm_pot_withdrawal_dates[date_index] = date
        
        # add the value to the dictionary
        self.groundwater_longterm_potential_withdrawal[date] = groundwater_longterm_pot_withdrawal
        
        # echo to screen
        message_str = str.join(' ', \
                                ('groundwater potential withdrawals updated for', \
                                message_str))
        logger.debug(message_str)
        
        # [ surface water ]
        #if self.pumping_capacity_flag['surfacewater']:
        # get surface water potential withdrawals for date
        # (units: m3/day)
        surfacewater_potential_withdrawal = self.surfacewater_potential_estimated_withdrawal
        
        # get the time step to update the groundwater potential withdrawal (monthly)
        # variable matches with groundwater_longterm_avail_dates
        date_index, matched_date, message_str = \
                    match_date_by_julian_number(date, \
                                                self.surfacewater_longterm_pot_withdrawal_dates)
        
        # remove the date from the dictionary and update it with the present value
        # set the value using the weight, if the long-term availability is not
        # defined, cover with the present value
        surfacewater_longterm_pot_withdrawal = self.surfacewater_longterm_potential_withdrawal.pop(matched_date)
        surfacewater_longterm_pot_withdrawal = \
              pcr.cover( self.surfacewater_update_weight      * surfacewater_potential_withdrawal + \
                        (1 - self.surfacewater_update_weight) * surfacewater_longterm_pot_withdrawal, \
                        surfacewater_potential_withdrawal)
        
        # reset the date
        self.surfacewater_longterm_pot_withdrawal_dates[date_index] = date
        
        # add the value to the dictionary
        self.surfacewater_longterm_potential_withdrawal[date] = surfacewater_longterm_pot_withdrawal
        
        # echo to screen
        message_str = str.join(' ', \
                                ('surface water potential withdrawals updated for', \
                                message_str))
        logger.debug(message_str)
        
        # returns None
        return None
    
    
    
    def get_final_conditions(self):

        '''
        get_final_conditions: function that returns a dictionary holding all states and fluxes
                             of the soil hydrology module that are necessary for a restart.
        
        output:
        ======
        state_info         : a dictionary with the key and the value
        '''
        
        # initialize states
        state_info = {}
        
        # iterate over the report name and attribute name
        for report_name, attr_name in self.report_state_info.items():
            
            # set the state
            state_info[report_name] = getattr(self, attr_name)
        
        # return results
        return state_info

# ///  end of the water management class ///
