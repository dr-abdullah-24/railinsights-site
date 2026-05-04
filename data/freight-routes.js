// Network Rail Intermodal Freight Terminals & Key Corridors
// Source: Network Rail freight map - intermodal sector (2021/2022)
window.FREIGHT_DATA = {

// ── TERMINALS ─────────────────────────────────────────────────────────────
// [id, name, type, status, lat, lon]
// type: SRFI=Strategic Rail Freight Interchange, RFI=Rail Freight Interchange
// status: A=Active, O=Operational(infrequent), U=Under development, P=Proposed
terminals: [
  [1,  "Inverness",                       "RFI",  "A", 57.4809,-4.2241],
  [2,  "Aberdeen",                         "RFI",  "A", 57.1437,-2.0981],
  [3,  "Grangemouth (Forth Ports)",        "RFI",  "A", 56.0167,-3.7167],
  [4,  "Grangemouth (DB Cargo)",           "RFI",  "A", 56.0183,-3.7150],
  [5,  "Grangemouth (Malcolm Group)",      "RFI",  "A", 56.0200,-3.7133],
  [6,  "Coatbridge",                       "RFI",  "A", 55.8600,-4.0283],
  [7,  "Mossend Eurocentral",              "SRFI", "A", 55.8283,-3.9817],
  [8,  "Mossend International",            "SRFI", "A", 55.8250,-3.9833],
  [9,  "Elderslie",                        "RFI",  "A", 55.8467,-4.4700],
  [10, "Workington",                       "RFI",  "O", 54.6440,-3.5444],
  [11, "Middlesbrough",                    "RFI",  "A", 54.5742,-1.2347],
  [12, "Teesport",                         "RFI",  "O", 54.5933,-1.1317],
  [13, "Seaforth (Liverpool)",             "RFI",  "A", 53.4500,-3.0200],
  [14, "Garston (Liverpool)",              "RFI",  "A", 53.3500,-2.8800],
  [15, "Widnes",                           "RFI",  "A", 53.3650,-2.7300],
  [16, "Parkside East",                    "SRFI", "P", 53.4550,-2.5500],
  [17, "Port Salford",                     "SRFI", "U", 53.4783,-2.3333],
  [18, "Trafford Park (DB Cargo)",         "RFI",  "A", 53.4617,-2.3217],
  [19, "Trafford Park (Freightliner)",     "RFI",  "A", 53.4600,-2.3200],
  [20, "Sheffield Masborough",             "RFI",  "A", 53.4100,-1.3617],
  [21, "Sheffield Tinsley",               "RFI",  "A", 53.4167,-1.3833],
  [22, "iPort Doncaster",                 "SRFI", "A", 53.5700,-0.9800],
  [23, "Doncaster Railport",              "RFI",  "A", 53.5250,-1.1333],
  [24, "Wakefield Europort",              "SRFI", "A", 53.6700,-1.4583],
  [25, "Leeds Stourton",                  "RFI",  "A", 53.7700,-1.5167],
  [26, "Konect",                          "RFI",  "U", 53.7833,-0.9833],
  [27, "Gascoigne Interchange",           "SRFI", "A", 53.7900,-1.2367],
  [28, "Hull King George Dock",           "RFI",  "O", 53.7450,-0.3017],
  [29, "Immingham",                       "RFI",  "O", 53.6267,-0.1867],
  [30, "Telford",                         "RFI",  "O", 52.6750,-2.4417],
  [31, "West Midlands Interchange",       "SRFI", "U", 52.6617,-2.2333],
  [32, "Cannock",                         "RFI",  "P", 52.6900,-2.0117],
  [33, "Burton on Trent",                 "RFI",  "O", 52.8150,-1.6350],
  [34, "East Midlands Intermodal Park",   "SRFI", "P", 52.8317,-1.3567],
  [35, "Castle Donington EMDC",           "RFI",  "O", 52.8283,-1.3583],
  [36, "East Midlands Gateway",           "SRFI", "A", 52.8267,-1.2667],
  [37, "Birmingham Lawley St",            "RFI",  "A", 52.4883,-1.8617],
  [38, "Hams Hall",                       "SRFI", "A", 52.5233,-1.7317],
  [39, "Birch Coppice",                   "SRFI", "A", 52.5850,-1.6667],
  [40, "Hinckley",                        "SRFI", "P", 52.5417,-1.3733],
  [41, "DIRFT I (Malcolm Group)",         "SRFI", "A", 52.3383,-1.1717],
  [42, "DIRFT II (Sainsbury's)",          "SRFI", "A", 52.3333,-1.1583],
  [43, "DIRFT II (Tesco)",               "SRFI", "A", 52.3350,-1.1650],
  [44, "DIRFT III",                       "SRFI", "U", 52.3300,-1.1533],
  [45, "Rail Central",                    "SRFI", "P", 52.1750,-1.0283],
  [46, "Northampton Gateway",             "SRFI", "U", 52.2167,-0.9500],
  [47, "Ely",                             "RFI",  "O", 52.3983, 0.2617],
  [48, "Barry Docks",                     "RFI",  "A", 51.3950,-3.2600],
  [49, "Wentloog (Cardiff)",              "RFI",  "A", 51.5017,-3.1100],
  [50, "Portbury",                        "RFI",  "O", 51.4950,-2.7433],
  [51, "Avonmouth",                       "RFI",  "A", 51.5100,-2.7250],
  [52, "Bristol South Liberty Lane",      "RFI",  "O", 51.4283,-2.5417],
  [53, "Swindon Keypoint",               "RFI",  "O", 51.5483,-1.7683],
  [54, "Sundon Quarry (Luton)",          "RFI",  "P", 51.9100,-0.4083],
  [55, "Radlett",                         "SRFI", "U", 51.7167,-0.3150],
  [56, "Willesden Euroterminal",          "RFI",  "O", 51.5317,-0.2483],
  [57, "Barking",                         "RFI",  "A", 51.5400, 0.0750],
  [58, "Dagenham (Ford)",                 "RFI",  "A", 51.5317, 0.1450],
  [59, "Purfleet",                        "RFI",  "A", 51.4767, 0.2433],
  [60, "Tilbury 1a",                      "RFI",  "A", 51.4683, 0.3550],
  [61, "Tilbury 1b",                      "RFI",  "A", 51.4700, 0.3567],
  [62, "Tilbury 2",                       "RFI",  "A", 51.4717, 0.3600],
  [63, "London Gateway",                  "RFI",  "A", 51.5067, 0.4667],
  [64, "Thames Enterprise Park",          "RFI",  "U", 51.5100, 0.5200],
  [65, "Thamesport",                      "RFI",  "O", 51.4217, 0.6033],
  [66, "Harwich",                         "RFI",  "A", 51.9467, 1.2833],
  [67, "Ipswich Griffin Wharf",           "RFI",  "A", 52.0500, 1.1533],
  [68, "Felixstowe",                      "RFI",  "A", 51.9633, 1.3517],
  [69, "Felixstowe North",               "RFI",  "A", 51.9650, 1.3500],
  [70, "Felixstowe South",               "RFI",  "A", 51.9617, 1.3533],
  [71, "Southampton Maritime",            "RFI",  "A", 50.9017,-1.4050],
  [72, "Southampton Millbrook",           "RFI",  "A", 50.9000,-1.4267],
  [73, "Southampton 107-108 Berth",       "RFI",  "A", 50.8983,-1.4183],
  [74, "Fratton (Portsmouth)",            "RFI",  "O", 50.8000,-1.0717]
],

// ── KEY FREIGHT CORRIDORS ─────────────────────────────────────────────────
// gauge: W12=highest, W10=standard, W9=Channel Tunnel route
// Each corridor: { name, gauge, colour, coords:[[lat,lon],...] }
corridors: [
  {
    name:"Felixstowe–Nuneaton Corridor",
    gauge:"W12", colour:"#00796b",
    info:"25kV AC electrified (Ipswich–Peterborough). W12 cleared. Key container route from Port of Felixstowe.",
    coords:[
      [51.963,1.352],[51.966,1.302],[52.050,1.153],[52.143,1.065],
      [52.187,0.977],[52.283,0.541],[52.398,0.262],
      [52.571,0.241],[52.574,0.237],
      [52.574,0.243],[52.648,-0.128],[52.577,-0.241],
      [52.570,-0.247],[52.574,-0.481],[52.666,-0.624],
      [52.777,-0.741],[52.956,-0.912],[53.032,-0.947],
      [53.085,-1.017],[53.160,-1.123],
      [52.960,-1.150],[52.660,-1.200],[52.400,-1.100],
      [52.338,-1.168]
    ]
  },
  {
    name:"WCML – West Coast Main Line (Freight)",
    gauge:"W10", colour:"#1565c0",
    info:"25kV AC OLE electrified. W10 gauge. UK's busiest mixed-use freight corridor. Max 75mph freight.",
    coords:[
      [51.528,-0.133],[51.534,-0.169],[51.541,-0.292],
      [51.619,-0.424],[51.702,-0.425],
      [51.855,-0.410],[52.081,-0.552],[52.164,-0.754],
      [52.244,-0.884],[52.310,-1.097],[52.338,-1.168],
      [52.380,-1.508],[52.427,-1.645],[52.468,-1.742],
      [52.514,-1.943],[52.617,-1.960],[52.690,-2.050],
      [52.790,-2.123],[53.020,-2.195],[53.078,-2.203],
      [53.260,-2.452],[53.373,-2.698],[53.430,-2.961],
      [53.750,-2.717],[53.870,-2.608],[54.010,-2.539],
      [54.119,-2.614],[54.205,-2.742],[54.328,-2.923],
      [54.550,-2.939],[54.688,-2.880],[54.900,-3.015],
      [55.090,-3.345],[55.350,-3.563],[55.573,-3.710],
      [55.752,-3.946],[55.866,-4.110],[55.861,-4.242]
    ]
  },
  {
    name:"ECML – East Coast Main Line (Freight)",
    gauge:"W10", colour:"#1976d2",
    info:"25kV AC OLE electrified. W10 gauge. Links London to Yorkshire, Humber, Tees and Scotland.",
    coords:[
      [51.533,-0.123],[51.576,-0.098],[51.726,-0.105],
      [51.887,-0.195],[52.105,-0.216],[52.340, 0.007],
      [52.574,-0.241],[52.806,-0.372],[52.960,-0.481],
      [53.062,-0.772],[53.103,-0.937],[53.160,-1.123],
      [53.302,-1.100],[53.518,-1.091],[53.521,-1.133],
      [53.697,-1.340],[53.775,-1.549],[53.834,-1.545],
      [53.958,-1.096],[54.001,-1.086],[54.188,-1.262],
      [54.404,-1.547],[54.545,-1.552],[54.671,-1.376],
      [54.769,-1.572],[54.841,-1.564],[54.975,-1.613],
      [55.052,-1.707],[55.184,-1.756],[55.389,-1.819],
      [55.502,-1.867],[55.627,-2.074],[55.747,-2.402],
      [55.870,-2.660],[55.938,-3.004],[55.957,-3.186],
      [56.000,-3.347],[56.038,-3.432]
    ]
  },
  {
    name:"Southampton–Midlands Corridor",
    gauge:"W10", colour:"#6a1599",
    info:"W10 gauge. Connects Port of Southampton to Midlands via Reading and Didcot. Key automotive and container route.",
    coords:[
      [50.901,-1.405],[50.935,-1.372],[50.985,-1.305],
      [51.076,-1.301],[51.136,-1.234],[51.174,-1.159],
      [51.263,-1.087],[51.385,-1.006],[51.458,-1.011],
      [51.457,-0.974],[51.465,-0.954],[51.505,-1.004],
      [51.575,-1.002],[51.604,-0.988],[51.665,-1.100],
      [51.716,-1.234],[51.750,-1.241],[51.762,-1.302],
      [51.820,-1.404],[51.894,-1.451],[52.064,-1.614],
      [52.204,-1.692],[52.324,-1.726],[52.427,-1.645],
      [52.468,-1.742],[52.514,-1.943]
    ]
  },
  {
    name:"Felixstowe–London via Ipswich",
    gauge:"W10", colour:"#2e7d32",
    info:"W10 gauge. Route from Felixstowe to London via Ipswich, Colchester, Chelmsford.",
    coords:[
      [51.963,1.352],[52.050,1.153],[52.060,1.143],
      [51.998,0.953],[51.930,0.873],[51.760,0.700],
      [51.731,0.473],[51.660,0.410],[51.540,0.075],
      [51.534,-0.017],[51.528,-0.123]
    ]
  },
  {
    name:"Trans-Pennine Freight Corridor",
    gauge:"W10", colour:"#e65100",
    info:"W10 gauge. Links Liverpool, Manchester, Sheffield and Humber ports across the Pennines.",
    coords:[
      [53.430,-2.961],[53.480,-2.897],[53.479,-2.334],
      [53.462,-2.322],[53.476,-2.218],[53.487,-2.083],
      [53.496,-1.997],[53.498,-1.916],[53.497,-1.770],
      [53.497,-1.616],[53.497,-1.549],[53.491,-1.388],
      [53.410,-1.362],[53.360,-1.182],[53.302,-1.100],
      [53.160,-1.123],[53.085,-1.017],[53.051,-0.861],
      [53.032,-0.803],[52.993,-0.608],[53.051,-0.329],
      [53.162,-0.185],[53.279,-0.112],[53.404,-0.031],
      [53.627,-0.187]
    ]
  },
  {
    name:"Scotland – ECML Extension",
    gauge:"W10", colour:"#00838f",
    info:"25kV AC OLE. Route linking Edinburgh to Glasgow, Coatbridge/Mossend, Grangemouth and beyond.",
    coords:[
      [56.000,-3.347],[55.957,-3.186],[55.938,-3.004],
      [55.947,-3.190],[55.950,-3.348],[55.957,-3.455],
      [55.950,-3.650],[55.866,-3.940],[55.861,-4.020],
      [55.866,-4.110],[55.861,-4.242]
    ]
  },
  {
    name:"Scotland – Inverness Corridor",
    gauge:"W10", colour:"#00695c",
    info:"Key Scottish corridor linking Aberdeen, Inverness and the central belt.",
    coords:[
      [55.950,-3.650],[56.038,-3.432],[56.119,-3.385],
      [56.404,-3.433],[56.468,-3.011],[56.590,-2.873],
      [56.710,-2.641],[56.902,-2.458],[57.143,-2.098],
      [57.242,-2.158],[57.323,-2.363],[57.384,-2.682],
      [57.481,-4.224]
    ]
  },
  {
    name:"Bristol/Cardiff–Midlands Corridor",
    gauge:"W10", colour:"#ad1457",
    info:"W10 gauge. Connects South Wales and Bristol ports to Midlands via Gloucester and Cheltenham.",
    coords:[
      [51.395,-3.260],[51.465,-3.120],[51.502,-3.110],
      [51.510,-2.725],[51.495,-2.743],[51.470,-2.617],
      [51.428,-2.542],[51.453,-2.580],[51.458,-2.604],
      [51.549,-1.768],[51.716,-1.962],[51.755,-1.987],
      [51.883,-2.073],[51.882,-2.073],[51.946,-2.078],
      [52.058,-2.097],[52.190,-2.217],[52.324,-1.726],
      [52.427,-1.645],[52.468,-1.742]
    ]
  },
  {
    name:"HS1 / Channel Tunnel – London",
    gauge:"W9", colour:"#c62828",
    info:"W9/UIC GB1 loading gauge. 25kV AC electrified. High Speed 1 line from Channel Tunnel to London St Pancras. Max 300km/h passenger, 120km/h freight.",
    coords:[
      [51.165, 1.123],[51.153, 1.057],[51.137, 0.939],
      [51.310, 0.654],[51.358, 0.540],[51.403, 0.440],
      [51.428, 0.305],[51.445, 0.248],[51.445, 0.122],
      [51.442,-0.017],[51.515, 0.013],[51.532,-0.023],[51.531,-0.124]
    ]
  },
  {
    name:"Humber Ports Corridor",
    gauge:"W10", colour:"#37474f",
    info:"Serves Immingham and Hull ports. Key bulk and container freight corridor.",
    coords:[
      [53.627,-0.187],[53.627,-0.186],[53.626,-0.187],
      [53.745,-0.302],[53.741,-0.301],
      [53.697,-0.430],[53.697,-1.340],[53.775,-1.549]
    ]
  }
],

// ── STATUS COLOURS ────────────────────────────────────────────────────────
statusColour: { A:"#0a7c4c", O:"#b45309", U:"#1565c0", P:"#6b7a8e" },
statusLabel:  { A:"Active", O:"Operational (infrequent)", U:"Under development", P:"Proposed" },
typeLabel:    { SRFI:"Strategic Rail Freight Interchange", RFI:"Rail Freight Interchange" },

// ── GAUGE LEGEND ─────────────────────────────────────────────────────────
gaugeLegend: [
  { gauge:"W12", colour:"#00796b", desc:"W12 – highest clearance, allows 9ft 6in Hi-Cube containers" },
  { gauge:"W10", colour:"#1565c0", desc:"W10 – standard UK freight gauge, 8ft 6in containers" },
  { gauge:"W9",  colour:"#c62828", desc:"W9 / UIC GB1 – Channel Tunnel route (HS1)" }
]

};
