var overlays = [];
var markers = [];
var map = null;
var currentOverlay = "all";

function addButtonFor(station, val) {
    var template = document.getElementById("station-template");
    var clone = template.cloneNode(true);
    clone.classList.remove("hidden");
    var button = clone.querySelector("button");
    button.innerHTML = station;
    button.addEventListener("click", selectOverlay.bind(undefined, val));
    document.getElementById("station-table-body").appendChild(clone);
}

function addBlankRow() {
    var template = document.getElementById("station-template");
    var clone = template.cloneNode(true);
    clone.classList.remove("hidden");
    var button = clone.querySelector("button");
    button.parentNode.removeChild(button);
    document.getElementById("station-table-body").appendChild(clone);
}

function initialize() {
    document.getElementById("date_start").innerHTML = first_position;
    document.getElementById("date_end").innerHTML = last_position;
    document.getElementById("num_pos").innerHTML = num_positions;

    var absbounds = null;

    addButtonFor("All coverage", "all");
    addBlankRow();
    addButtonFor("4+ station overlap", "4plus");
    addButtonFor("5+ station overlap", "5plus");
    addButtonFor("6+ station overlap", "6plus");
    addButtonFor("Below 18000ft", "below18000");
    addButtonFor("Below 10000ft", "below10000");
    addButtonFor("Min altitude seen", "byalt");
    addBlankRow();
    
    // Create map
    map = L.map('map-canvas');

    // create the OSM tile layer
    var osmUrl='https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png';
    var osmAttrib='Map data <a href="https://openstreetmap.org">OpenStreetMap</a> contributors';
    var osm = new L.TileLayer(osmUrl, {minZoom: 2, maxZoom: 20, attribution: osmAttrib});

    // default map start position
    map.setView(new L.LatLng(51.781, -1.5706), 6);

    // Add the OSM layer to the map
    map.addLayer(osm);
        
    var names = Object.keys(coverage).sort();
    for (var i = 0; i < names.length; ++i) {
        var k = names[i];
        var s = coverage[k];
        var bounds = new L.LatLngBounds(
            new L.LatLng(s.min_lat, s.min_lon),
            new L.LatLng(s.max_lat, s.max_lon));
                        
        overlays[k] = L.imageOverlay(s.image, bounds, { opacity : 1.0 });
        
        if (s.lat !== null) {
            // marker jitter is just to separate markers that would otherwise be overlapping
            markers[k] = L.marker([s.lat + Math.random()*0.02-0.01, s.lon + Math.random()*0.02-0.01]).addTo(map).bindTooltip(k);
            markers[k].on('click', selectOverlay.bind(undefined, k));
        }
                
        if (s.is_station) {
            addButtonFor(k, k);
        }
    }

    overlays['all'].addTo(map);
      
    var absbounds = new L.LatLngBounds(new L.LatLng(coverage['all'].min_lat, coverage['all'].min_lon), new L.LatLng(coverage['all'].max_lat, coverage['all'].max_lon));
       
    map.fitBounds(absbounds);
}

function selectOverlay(stationname) {
    map.removeLayer(overlays[currentOverlay]);

    if (currentOverlay === stationname) {
        stationname = "all";
    }
        
    overlays[stationname].addTo(map);
    currentOverlay = stationname;
}

