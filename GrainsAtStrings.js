/*
	chainedgsys.js
	Simulation of chained granular particles system
	
	Sparisoma Viridi | https://github.com/dudung/butiran
	Usman Sambiri | usmansucces@gmail.com 
	Inayatul Inayah | inayatul.inayah94@gmail.com
	Hazrati Ashel | hazratiashel18@gmail.com
	Zulhendra | zulhendraphetir@gmail.com
	Dellia Yulita | dellia.yulita@gmail.com
	Nofri ermasari | nofri.ermasari@gmail.com
	
	Include: <script src="chainedgsys.js"></script> in a HTML
	         file
	Execute: Refresh web browser viewing the HTML file
	
	20190412
	0537 Beg zuhause.
	0641 Con an der Uni.
	1604 Fin numberOf().
	20190413
	0805 Start again at Aula Barat.
	1049 Fix rounding of time t to tt.
	1126 Back to my cube and con.
	1132 Beg transform().
	1150 Beg visualize().
	1257 Add length of each spring.
	1600 Start at home.
	1604 Change fluid to water.
	1708 Still error.
	20190414
	0430 Fix problem without vertical forces with
	.toExponential(3) in dx and dy.
	1043 Con at smansa.
	1534 Fix explosion by introduction of spring damping.
	1544 Add authors for further development.
	
	References
	1 "How can I pad a value with leading zeros?"
	  url https://stackoverflow.com/a/20460414/9475509
		[20190413]
	2 "Round to at most 2 decimal places (only if necessary)"
	  url https://stackoverflow.com/a/12830454/9475509
		[20190413]
	3 Dynamic viscosity at 25 oC of water 8.9E-4
	  url https://www.engineersedge.com/physics
		/water__density_viscosity_specific_weight_13146.htm
		[20190413]
	4 Dynamic viscosity at 25 oC of air 18.6E-6
	  url https://www.engineersedge.com/physics
		/viscosity_of_air_dynamic_and_kinematic_14483.htm
		[20190413]
	5 At 20 °C and 101.325 kPa, dry air has a density of
	  1.2041 kg/m³
		url https://en.wikipedia.org/wiki/Density_of_air
		[20190414]
*/

// Define global variables for visual elements
var taIn, taOut, btRead, btStart, caOut;
var proc, Tproc, Nproc, iproc;

// Define coordinates
var xmin, ymin, xmax, ymax;
var XMIN, YMIN, XMAX, YMAX;

// Define physical parameters
var tbeg, tend, t, dt;
var g, eta, rho, T;
var x, y, vx, vy, m, D;
var k;
var Ngrains, Nspring;
var Us, Ug, K, E;

// Execute main function
main();


// Define main function
function main() {
	initElements();
	initParams();
	writeDefaultParams();
}


// Write default parameters to textarea
function writeDefaultParams() {
	// Write to Textarea
	tin("# Simulation\n");
	tin("TBEG 0\n");
	tin("TEND 10\n");
	tin("TSTEP 0.001\n");
	tin("\n");
	
	tin("# World coordinates\n");
	tin("XMIN -0.5\n");
	tin("YMIN -1.5\n");
	tin("XMAX 2.5\n");
	tin("YMAX 1.5\n");
	tin("\n");
	
	tin("# Environment\n");
	tin("GACC 9.81\n");     // m/s^2
	tin("ETAF 18.6E-6\n");  // Pa.s
	tin("RHOF 1.2041\n");   // kg/m^3
	tin("TENV 298\n");      // K
	tin("\n");
	
	tin("# Grains x y vx vy m D\n");
	tin("GRAIN00 0.0 0.0 0 0 0.6 0.1\n");
	tin("GRAIN01 0.2 0.0 0 0 0.6 0.1\n");
	tin("GRAIN02 0.4 0.0 0 0 0.6 0.1\n");
	tin("GRAIN03 0.6 0.0 0 0 0.6 0.1\n");
	tin("GRAIN04 0.8 0.0 0 0 0.6 0.1\n");
	tin("GRAIN05 1.0 0.0 0 0 0.6 0.1\n");
	tin("GRAIN06 1.2 0.0 0 0 0.6 0.1\n");
	tin("GRAIN07 1.4 0.0 0 0 0.6 0.1\n");
	tin("GRAIN08 1.6 0.0 0 0 0.6 0.1\n");
	tin("GRAIN09 1.8 0.0 0 0 0.6 0.1\n");
	tin("GRAIN10 2.0 0.0 0 0 0.6 0.1\n");
	tin("\n");
	
	tin("# Spring k l beg end\n");
	tin("SPRING0 2000 100 0.2 0 1\n");
	tin("SPRING1 2000 100 0.2 1 2\n");
	tin("SPRING2 2000 100 0.2 2 3\n");
	tin("SPRING3 2000 100 0.2 3 4\n");
	tin("SPRING4 2000 100 0.2 3 4\n");
	tin("SPRING5 2000 100 0.2 3 4\n");
	tin("SPRING6 2000 100 0.2 3 4\n");
	tin("SPRING7 2000 100 0.2 3 4\n");
	tin("SPRING8 2000 100 0.2 3 4\n");
	tin("SPRING9 2000 100 0.2 3 4\n");
	tin("\n");
}


// Initialize parameters
function initParams() {
	// Iteration parameters
	Tproc = 1;
	Nproc = 20;
	iproc = 0;
	
	// Canvas coordinates
	XMIN = 0;
	YMIN = caOut.height;
	XMAX = caOut.width;
	YMAX = 0;
}


// Initialize elements
function initElements() {
	taIn = document.createElement("textarea");
	taIn.style.width = "240px";
	taIn.style.height = "350px";
	taIn.style.overflowY = "scroll";
	taIn.style.float = "left";
	
	btStart = document.createElement("button");
	btStart.innerHTML = "Start";
	btStart.style.width = "50px";
	btStart.addEventListener("click", buttonClick);
	btStart.disabled = true;
	
	btRead = document.createElement("button");
	btRead.innerHTML = "Read";
	btRead.style.width = "50px";
	btRead.addEventListener("click", buttonClick);
	
	var divButton = document.createElement("div");
	divButton.style.width = "50px";
	divButton.style.float = "left";	
	
	taOut = document.createElement("textarea");
	taOut.style.width = "200px";
	taOut.style.height = "350px";
	taOut.style.overflowY = "scroll";
	taOut.style.float = "left";
	
	caOut = document.createElement("canvas");
	caOut.width = "354";
	caOut.height = "354";
	caOut.style.width = caOut.width;
	caOut.style.height = caOut.height;
	caOut.style.border = "1px solid #ccc";
	caOut.style.background = "#fafafa";
	
	document.body.append(taIn);
	document.body.append(divButton);
		divButton.append(btRead);
		divButton.append(btStart);
	document.body.append(taOut);
	document.body.append(caOut);
}


// Check number of same prefix keywords
function numberOf(keyword) {
	var val = 0;
	var lines = taIn.value.split("\n");
	for(var l = 0; l < lines.length; l++) {
		var words = lines[l].split(" ");
		if(keyword == words[0].substring(0, keyword.length)) {
			val++;
		}
	}
	return val;
}


// Read a parameter
function readParam(keyword) {
	var val;
	var lines = taIn.value.split("\n");
	for(var l = 0; l < lines.length; l++) {
		var words = lines[l].split(" ");
		if(keyword == words[0]) {
			if(words.length == 2) {
				val = parseFloat(words[1]);
			} else {
				val = [];
				for(var w = 1; w < words.length; w++) {
					var v = parseFloat(words[w]);
					val.push(v);
				}
			}
		}
	}
	return val;
}


// Read all parameters
function readAllParams() {
	//var xmin, ymin, xmax, ymax;
	xmin = readParam("XMIN"); tout("XMIN " + xmin + "\n");
	ymin = readParam("YMIN"); tout("YMIN " + ymin + "\n");
	xmax = readParam("XMAX"); tout("XMAX " + xmax + "\n");
	ymax = readParam("YMAX"); tout("YMAX " + ymax + "\n");
	
	//var tbeg, tend, t, dt;
	tbeg = readParam("TBEG"); tout("TBEG " + tbeg + "\n");
	tend = readParam("TEND"); tout("TEND " + tend + "\n");
	dt = readParam("TSTEP"); tout("TSTEP " + dt + "\n");
	
	Nproc = Math.ceil((tend - tbeg) / dt);
	
	//var g, eta, T;
	g = readParam("GACC"); tout("GACC " + g + "\n");
	eta = readParam("ETAF"); tout("ETAF " + eta + "\n");
	rho = readParam("RHOF"); tout("RHOF " + rho + "\n");
	T = readParam("TENV"); tout("TENV " + T + "\n");
	
	//var x, y, vx, vy, m, D;
	x = [];
	y = [];
	vx = [];
	vy = [];
	m = [];
	D = [];

	NGrains = numberOf("GRAIN");
	var digGrains = 1 + Math.floor(Math.log10(NGrains - 1));
	
	for(var i = 0; i < NGrains; i++) {
		var num = ("00000000" + i).slice(-digGrains);
		var keyword = "GRAIN" + num;
		var read = readParam(keyword);
		x.push(read[0]);
		y.push(read[1]);
		vx.push(read[2]);
		vy.push(read[3]);
		m.push(read[4]);
		D.push(read[5]);
	}
	
	//var k;
	k = [];
	
	NSprings = numberOf("SPRING");
	var digSprings = 1 + Math.floor(Math.log10(NSprings - 1));
	
	for(var i = 0; i < NSprings; i++) {
		var num = ("00000000" + i).slice(-digSprings);
		var keyword = "SPRING" + num;
		var read = readParam(keyword);
		k.push(read);
	}
}


// Do something when button is clicked
function buttonClick() {
	var t = event.target;
	if(t.innerHTML == "Start") {
		proc = setInterval(simulate, Tproc);
		t.innerHTML = "Stop";
		btRead.disabled = true;
	} else if(t.innerHTML == "Stop"){
		clearInterval(proc);
		t.innerHTML = "Start"
		btRead.disabled = false;
	} else if(t.innerHTML == "Read") {
		btStart.disabled = false;
		initParams();
		readAllParams();		
	}
}


// Perform calculation
function calculate() {
	// Calculate force for all grains except the first (0)
	// and the last (N - 1)
	
	// Define temporary value
	var xnew = x.slice();
	var ynew = y.slice();
	var vxnew = vx.slice();
	var vynew = vy.slice();
	
	// Calculate current total energy
	var Ecur = 0;
	var Ugcur = 0;
	var Uscur = 0;
	var Kcur = 0;
	for(var i = 1; i < NGrains-1; i++) {
		Ugcur += m[i] * g * y[i];
		
		for(var j = i-1; j < i+2; j += 2) {
			var dx = x[i] - x[j];
			var dy = y[i] - y[j];
			var dr = Math.sqrt(dx * dx + dy * dy);
			
			Uscur += 0.5 * k[i][0] * (dr - k[i][2]) * (dr - k[i][2]);
		}
		
		var vx2 = vx[i] * vx[i];
		var vy2 = vy[i] * vy[i];
		var v = Math.sqrt(vx2 + vy2);
		
		Kcur = 0.5 * m[i] * v * v;
		
		Ecur += (Ugcur + Uscur + Kcur);
	}
	
	// Calculate motion only using conservative forces
	for(var i = 1; i < NGrains-1; i++) {
		// Set variables for total force
		var Fx = 0;
		var Fy = 0;
		
		// Calculate gravitational force
		var Fgx = 0;
		var Fgy = m[i] * -g;
		
		Fx += Fgx;
		Fy += Fgy;
		
		// Calculate spring force from neighbour
		for(var j = i-1; j < i+2; j += 2) {
			var dx = (x[i] - x[j]).toExponential(10);
			//var dx = (x[i] - x[j]);
			var dy = (y[i] - y[j]).toExponential(10);
			//var dy = (y[i] - y[j]);
			var dr = Math.sqrt(dx * dx + dy * dy);
			//dx = Math.sqrt(dr * dr - dy * dy);
			//dy = Math.sqrt(dr * dr - dx * dx);
			//dr = Math.sqrt(dx * dx + dy * dy);
			
			var Fs = -k[i][0] * (dr - k[i][2]);
			var Fsx = (dr > 0) ? Fs * (dx / dr) : 0;
			var Fsy = (dr > 0) ? Fs * (dy / dr) : 0;
			//var dd2 = dr*dr - dx*dx - dy*dy;
			
			Fx += Fsx;
			Fy += Fsy;
			
			var dvx = vx[i] - vx[j];
			var dvy = vy[i] - vy[j];
			var dv = Math.sqrt(dvx * dvx + dvy * dvy);
			var Fsv = -k[i][1] * dv;
			var Fsvx = (dv > 0) ? Fsv * (dvx / dv) : 0;
			var Fsvy = (dv > 0) ? Fsv * (dvy / dv) : 0;
			
			Fx += Fsvx;
			Fy += Fsvy;
		}
		
		// Calculate buoyant force
		var R3 = (0.5 * D[i]) * (0.5 * D[i]) * (0.5 * D[i]);
		var V = (4.0 / 3.0) * Math.PI * R3;
		var Fb = rho * g * V;
		var Fbx = 0;
		var Fby = Fb;
		
		Fx += Fbx;
		Fy += Fby;
		

		// Calculate drag force
		var vx2 = vx[i] * vx[i];
		var vy2 = vy[i] * vy[i];
		var v = Math.sqrt(vx2 + vy2);
		var Fd = -3 * Math.PI * eta * D[i] * v;
		var Fdx = (v > 0) ? Fd * (vx[i] / v) : 0;
		var Fdy = (v > 0) ? Fd * (vy[i] / v) : 0;

		Fx += Fdx;
		Fy += Fdy;

		// Apply Newton's second law of motion
		var ax = Fx / m[i];
		var ay = Fy / m[i];
		
		// Integrate using Euler algoritm to get new velocity
		vxnew[i] = vx[i] + ax * dt;
		vynew[i] = vy[i] + ay * dt;
		
		// Integrate using Euler algoritm to get new position
		xnew[i] = x[i] + vx[i] * dt;
		ynew[i] = y[i] + vy[i] * dt;
	}
	
	// Calculate new total energy
	var Enew = 0;
	var Usnew = 0;
	var Ugnew = 0;
	var Knew = 0;
	for(var i = 1; i < NGrains-1; i++) {
		Ugnew += m[i] * g * ynew[i];
		
		for(var j = i-1; j < i+2; j += 2) {
			var dx = xnew[i] - xnew[j];
			var dy = ynew[i] - ynew[j];
			var dr = Math.sqrt(dx * dx + dy * dy);
			
			Usnew += 0.5 * k[i][0] * (dr - k[i][2]) * (dr - k[i][2]);
		}
		
		var vx2 = vxnew[i] * vxnew[i];
		var vy2 = vynew[i] * vynew[i];
		var v = Math.sqrt(vx2 + vy2);
		
		Knew = 0.5 * m[i] * v * v;
		
		Enew += (Ugnew + Usnew + Knew);
	}
	
	
	// Calculate increasing of kinetic energy
	var dE = Enew - Ecur;
	var dU = (Ugnew - Ugcur) + (Usnew - Uscur)
	var dK = dE - dU;
	
	K = Knew;
	Us = Usnew;
	Ug = Ugnew;
	E = Enew;
	
	/*
	Ecur = Ukur + Kkur
	Enew = Unew + Knew
	dE = Enew - Ecur;
	Knew = Knew - dK = Knew - dE
	     = Knew - ((Unew + Knew) - (Ucur + Kcur))
	     = Knew - Unew - Knew + Ucur + Kcur
			 = Kcur - (Unew - Ukur)
			 = Kcur - dU
	Kcor = Knew - dK
	     = Knew - (N - 2) dK / (N - 2)
	dKi = dK / (N - 2)
	Kicor = Ki - dKi
	vicor^2 = vi^2 - 2 dKi / mi
	vicor = (vi^2 - 2 dKi / mi)^0.5
	vicorx = vicor (vix/vi)
	vicory = vicor (viy/vi)
	*/
	
	/**/
	var dKi = dK / (NGrains - 2);
	for(var i = 1; i < NGrains-1; i++) {
		var vx2 = vxnew[i] * vxnew[i];
		var vy2 = vynew[i] * vynew[i];
		var v = Math.sqrt(vx2 + vy2);
		
		var vcor = Math.sqrt(v * v - 2 * dKi / m[i]);
		var vcorx = vcor * (vxnew[i] / v);
		var vcory = vcor * (vynew[i] / v);
		
		//vxnew[i] = vcorx;
		//vynew[i] = vcory;
	}
	/**/
	
	// Update all
	x = xnew.slice();
	y = ynew.slice();
	vx = vxnew.slice();
	vy = vynew.slice();	
}


// Perform visualisation
function visualize() {
	var cx = caOut.getContext("2d");
	cx.clearRect(0, 0, caOut.width, caOut.height);
	
	var N = NGrains;
	
	cx.strokeStyle = "#f00";
	cx.beginPath();
	for(var i = 0; i < N; i++) {
		var r = transform(x[i], y[i]);
		var X = r.x;
		var Y = r.y;
		if(i == 0) {
			cx.moveTo(X, Y);
		} else {
			cx.lineTo(X, Y);
		}
	}
	cx.stroke();
	cx.closePath();
	
	cx.fillStyle = "#ccf";
	cx.strokeStyle = "#00f";
	for(var i = 0; i < N; i++) {
		var r = transform(x[i], y[i]);
		var X = r.x;
		var Y = r.y;
		var r_dr = transform(x[i] - 0.5 * D[i], y[i]);
		var R = r.x - r_dr.x;
		cx.beginPath();
		cx.arc(X, Y, R, 0, 2 * Math.PI);
		cx.fill();
		cx.stroke();
		cx.closePath();
	}
}


// Perform calculation
function simulate() {
	if(iproc <= Nproc) {
		
		// Perform visualisation from previous iteration
		visualize();
		
		// Perform calculation for next iteration
		calculate();
		
		// Show time t
		t = iproc * dt;
		var tt = t.toFixed(-Math.floor(Math.log10(dt)));
		tout(
			iproc + " "
			+ tt + " "
			//+ Ug.toExponential(3) + " "
			//+ Us.toExponential(3) + " "
			//+ K.toExponential(3) + " "
			+ E.toExponential(5) + " "
			+ "\n"
		);
		
		// Advance iteration
		iproc++;
	} else {
		clearInterval(proc);
		btStart.innerHTML = "Start";
		btStart.disabled = true;
		btRead.disabled = false;
	}
}


// Transform coodinates from world to canvas
function transform(x, y) {
	var X = (x - xmin) / (xmax - xmin) * (XMAX - XMIN);
		X += XMIN;
	
	var Y = (y - ymin) / (ymax - ymin) * (YMAX - YMIN);
		Y += YMIN;
	
	return {x: X, y: Y};
}


// Add text output textarea
function tout(string) {
	taOut.value += string;
	taOut.scrollTop = taOut.scrollHeight;
}


// Add text input textarea
function tin(string) {
	taIn.value += string;
	taIn.scrollTop = taIn.scrollHeight;
}