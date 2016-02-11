function arrayUnique(arr) 
{
    return arr.reduce(function(p, c) 
    {
        if (p.indexOf(c) < 0) p.push(c);
        return p;
    }, []);
}

function setBoldExclusively(elem)
{
	$('.bold').removeClass('bold');
	elem.addClass('bold');
}

function shadeColor(color, percent) 
{   
    var f=parseInt(color.slice(1),16),t=percent<0?0:255,p=percent<0?percent*-1:percent,R=f>>16,G=f>>8&0x00FF,B=f&0x0000FF;
    return "#"+(0x1000000+(Math.round((t-R)*p)+R)*0x10000+(Math.round((t-G)*p)+G)*0x100+(Math.round((t-B)*p)+B)).toString(16).slice(1);
}

function blendColors(c0, c1, p) 
{
    var f=parseInt(c0.slice(1),16),t=parseInt(c1.slice(1),16),R1=f>>16,G1=f>>8&0x00FF,B1=f&0x0000FF,R2=t>>16,G2=t>>8&0x00FF,B2=t&0x0000FF;
    return "#"+(0x1000000+(Math.round((R2-R1)*p)+R1)*0x10000+(Math.round((G2-G1)*p)+G1)*0x100+(Math.round((B2-B1)*p)+B1)).toString(16).slice(1);
}

function onlyUnique(value, index, self) 
{ 
    return self.indexOf(value) === index;
}

function tabulateNumClusters(bootstraps)
{
	var result = new Array();
	result["connComponents"] = new Array();
	result["dipMeans"] = new Array();

	for (var i in bootstraps)
	{
		var bs = bootstraps[i];

		if (!result["connComponents"][bs.resConnComponents.numClusters])
		{			
			result["connComponents"][bs.resConnComponents.numClusters] = 0;
		}

		if (!result["dipMeans"][bs.resDipMeans.numClusters])
		{			
			result["dipMeans"][bs.resDipMeans.numClusters] = 0;
		}

		result["connComponents"][bs.resConnComponents.numClusters]++;
		result["dipMeans"][bs.resDipMeans.numClusters]++;
	}

	for (var i in result["connComponents"])
	{
		result["connComponents"][i] /= bootstraps.length;	
	}

	for (var i in result["dipMeans"])
	{
		result["dipMeans"][i] /= bootstraps.length;	
	}

	return result;
}

function calculateStats(results)
{
	stats = new Array();

	for (var i in results)
	{
		var res = results[i];
		
		var key = res.fasta;
		stats[key] = {};

		var tnc = tabulateNumClusters(res.bootstraps);
		stats[key].connComponents = tnc["connComponents"];
		stats[key].dipMeans = tnc["dipMeans"];

		stats[key].kraken = {};
		var krakenNumUnknown = 0;
		$.each(res.krakenLabels, function (idx, val) 
		{
			if (val === "unknown")
			{
				krakenNumUnknown++;
			}
		});
		stats[key].kraken.numUnknown = krakenNumUnknown;
		var krakenUnique = res.krakenLabels.filter(onlyUnique);
		stats[key].kraken.numSpecies = krakenUnique.length;
		if(krakenNumUnknown > 0)
		{
			stats[key].kraken.numSpecies--;
		}
	}

	return stats;
}

function cellBarChart(containerCell, confidences, maxK)
{
	// convert confidences into data array
	var data = new Array(maxK);
	for (var i = 1; i <= maxK; i++) 
	{
		data[i-1] = {};
		data[i-1].numClusters = i;
		if (i in confidences)
		{
			data[i-1].confidence = confidences[i];
		} else
		{
			data[i-1].confidence = 0;
		}
	}

	var margin = {top: 10, right: 20, bottom: 25, left: 20},
	width = 80 - margin.left - margin.right,
	height = 60 - margin.top - margin.bottom;

	var x = d3.scale.ordinal().rangeRoundBands([0, width], .1);
	var y = d3.scale.linear().range([height, 0]);

	var xAxis = d3.svg.axis().scale(x).orient("bottom");
	var yAxis = d3.svg.axis().scale(y).orient("left").ticks(1);

	var svg = d3.select($(containerCell).get(0)).append("svg")
		.attr("width", width + margin.left + margin.right)
		.attr("height", height + margin.top + margin.bottom)
		.append("g")
		.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	x.domain(data.map(function(d) { return d.numClusters; }));
	y.domain([0, 1]);


  	svg.append("g")
		.attr("class", "x axis")
		.attr("transform", "translate(0," + height + ")")
		.call(xAxis);

  	svg.append("g")
		.attr("class", "y axis")
		.call(yAxis);

    tip = d3.tip().attr('class', 'd3-tip').offset([-10, 0]).html(function(d, i) { return 'p(' + data[i].numClusters + ') = ' + data[i].confidence; });
    svg.call(tip);

	svg.selectAll(".bar")
		.data(data)
		.enter().append("rect")
		.attr("class", "bar")
		.attr("x", function(d) { return x(d.numClusters); })
		.attr("width", x.rangeBand())
		.attr("y", function(d) { return y(d.confidence); })
		.attr("height", function(d) { return height - y(d.confidence); })
		.on('mouseover', tip.show)
		.on('mouseout', tip.hide);

}

function buildConfidenceTable(results)
{
	var stats = calculateStats(results);

	var maxK = 0;
	for (var i in stats)
	{
		var keys =      Object.keys(stats[i].connComponents).map(Number)
				.concat(Object.keys(stats[i].dipMeans).map(Number));
				// .concat(stats[i].kraken.numSpecies);
		var maxKey = Math.max.apply(Math, keys);		
		maxK = Math.max(maxK, maxKey);
	}

	for (var i in stats)
	{
		// connected components

		var s = stats[i].connComponents;
		var cc = '';
		for (var k in s)
		{
			cc = cc + "p(k = " + k + ") = " + s[k] + "<br/>";
		}
		var cleanVal = stats[i].connComponents[1] ? stats[i].connComponents[1] : 0;
		var ccColor = blendColors('#f46d43', '#abdda4', cleanVal);

		// dip means
		
		var s = stats[i].dipMeans;
		var dip = '';
		for (var k in s)
		{
			dip = dip + "p(k = " + k + ") = " + s[k] + "<br/>";
		}
		var cleanVal = stats[i].dipMeans[1] ? stats[i].dipMeans[1] : 0;
		var dipColor = blendColors('#f46d43', '#abdda4', cleanVal);

		// kraken

		var kraken = stats[i].kraken.numSpecies + ' species';
		var krakenColor = stats[i].kraken.numSpecies > 1 ? '#f46d43' : '#abdda4';		

		$('#confidences').append('<tr><td class="dataConf">' + i + '</td><td style="background-color:' + ccColor + ';" class="ccConf"><div class="chart"></div></td><td style="background-color:' + dipColor + ';" class="dipConf"><div class="chart"></div></td><td style="background-color:' + krakenColor + ';" class="kraken">' + kraken + '</td></tr>');

		cellBarChart($('#confidences tr:last td.ccConf div.chart'), stats[i].connComponents, maxK);
		cellBarChart($('#confidences tr:last td.dipConf div.chart'), stats[i].dipMeans, maxK);
	}	
}

function dimMax(data, dim)
{
	return data.reduce(function(max, arr) 
	{
	    return max >= arr[dim] ? max : arr[dim];
	}, -Infinity);
}

function dimMin(data, dim)
{
	return data.reduce(function(min, arr) 
	{
	    return min < arr[dim] ? min : arr[dim];
	}, Infinity);
}

function showData(dataMat, labels, tooltips, width, height, padding)
{	
	var svg = d3.select("#scatter").html("").append("svg").attr("width", width).attr("height", height);

	svg.append("rect")
	   .attr("x", 0)
	   .attr("y", 0)
	   .attr("width", width)
	   .attr("height", height)
	   .style("fill", "none")
  	   .style("stroke", "lightgrey")
       .style("stroke-width", "1");	   

	minX = dimMin(dataMat, 0);
	maxX = dimMax(dataMat, 0);
	minY = dimMin(dataMat, 1);
	maxY = dimMax(dataMat, 1);

	var xScale = d3.scale.linear()
			.domain([minX, maxX])
			.range([padding/2, width-padding/2]);
	var yScale = d3.scale.linear()
             .domain([minY, maxY])
             .range([padding/2, height-padding/2]);

    tip = d3.tip().attr('class', 'd3-tip').offset([-10, 0]).html(function(d, i) { return tooltips[i]; });
    svg.call(tip);

	svg.selectAll("circle")
			.data(dataMat)
			.enter()
			.append("circle")
			.attr("cx", function(d) { return xScale(d[0]); })
       		.attr("cy", function(d) { return yScale(d[1]); })
       		.attr("r", function(d) {return 2; })
       		.style("fill", function(d, i) {return colors[labels[i] % colors.length]; })
			.on('mouseover', tip.show)
			.on('mouseout', tip.hide);	
}

function showVisualization()
{
	var x = results[selectedFasta];

	var clustAnaResult;
	if (selectedData === 'oneshot')
	{
		clustAnaResult = x.oneshot;
	} else
	{
		clustAnaResult = x.bootstraps[parseInt(selectedData)];
	}

	var dataMat = clustAnaResult[selectedReduction];

	var labels;
	if (selectedLabels === 'fasta')
	{
		// labels = numericLabels(x.fastaLabels);
		labels = Array.apply(null, Array(x.fastaLabels.length)).map(Number.prototype.valueOf, -1); // no labels / black color
	} else if(selectedLabels === 'cc')
	{
		labels = clustAnaResult.resConnComponents.labels;
	} else if(selectedLabels === 'dip')
	{
		labels = clustAnaResult.resDipMeans.labels;
	} else if(selectedLabels === 'kraken')
	{
		labels = numericLabels(x.krakenLabels);
	}

	var tooltips = results[selectedFasta].fastaLabels;
	if(selectedLabels === 'kraken')
	{
		tooltips = x.krakenLabels;
	}	

	showData(dataMat, labels, tooltips, width, height, padding);
	updateExport(labels);
}

function updateBootStrapsSelect()
{
	$('#bootstraps').val('oneshot');
	selectedData = 'oneshot';

	if (selectedLabels === 'kraken')
	{
		$('#bootstraps').prop("disabled", true);
		return;
	}

	$('#bootstraps').prop("disabled", false);

	var n = results[selectedFasta].bootstraps.length;

	$('#bootstraps').html('');
	$('#bootstraps').append('<option value="oneshot">One shot</option>');
	for (var i = 0; i < n; i++) 
	{
		$('#bootstraps').append($('<option>', {
		    value: i,
		    text: 'Bootstrap ' + (i+1)
		}));		
	}

}

function updateExport(labels)
{
	if (selectedLabels === 'fasta' || selectedData !== 'oneshot')
	{
		$('#export').hide();
		return;
	}

	var uniqueLabels = arrayUnique(labels).sort();

	$('#exportColors').html('');
	for (var i in uniqueLabels)
	{
		var lbl = uniqueLabels[i];
		var color = colors[lbl % colors.length];
		var href = 'export/' + results[selectedFasta].id + "-" + selectedLabels + "-" + lbl + ".fasta";
		$('#exportColors').append('<span><a style="color: ' + color + ';" href="' + href + '">&#x25CF;</a></span>&nbsp;');
	}

	$('#export').show();
}

function numericLabels(labels)
{
	var mp = new Array();
	k = 0;
	for (var i in labels) 
	{
		var key = labels[i];
		if (!(key in mp))
		{
			mp[key] = k;
			k++;
		}
	}
	var result = Array(labels.length);
	for (var i in labels) 
	{
		result[i] = mp[labels[i]];
	}
	return result;
}

function bootstrapLabels(labelsOneshot, bootstrapIndices)
{
	result = Array(bootstrapIndices.length);
	for (var i in bootstrapIndices)
	{
		result[i] = labelsOneshot[bootstrapIndices[i]];
	}	
	return result;
}
