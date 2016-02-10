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

function buildConfidenceTable(results)
{
	var stats = calculateStats(results);

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

		var kraken = 'Unique species: ' + stats[i].kraken.numSpecies + '<br/>Unknown: ' + stats[i].kraken.numUnknown;
		var krakenColor = stats[i].kraken.numSpecies > 1 ? '#f46d43' : '#abdda4';		

		$('#confidences').append('<tr><td class="dataConf">' + i + '</td><td class="ccConf" style="background-color:' + ccColor + '">' + cc + '</td><td class="dipConf" style="background-color:' + dipColor + '">' + dip + '</td><td class="kraken" style="background-color:' + krakenColor + '">' + kraken + '</td></tr>');
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
	var colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"];

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
}

function updateBootStrapsSelect()
{
	$('#bootstraps').val('oneshot');

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
