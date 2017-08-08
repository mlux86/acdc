var isFileSaverSupported = false;
try 
{
    var isFileSaverSupported = !!new Blob;
} catch (e) 
{
}
if(!isFileSaverSupported)
{
	alert('HTML5 Blobs are not supported and fasta export will not work. Please upgrade your browser.');
}

function stickyScatter()
{
    $('#scatterContainer').addClass('stick');
	var top = Math.max($(window).scrollTop() + 10, $('#confidences').offset().top);
    var left = $('#confidences').offset().left + $('#confidences').width() + 25;
	$('#scatterContainer').offset({ left: left, top: top });
}

function arrayUnique(arr)
{
    return arr.reduce(function(p, c)
    {
        if (p.indexOf(c) < 0) p.push(c);
        return p;
    }, []);
}

function arraySortGetIndices(toSort) 
{
	for (var i = 0; i < toSort.length; i++) 
	{
		toSort[i] = [toSort[i], i];
	}
	toSort.sort(function(left, right) 
	{
		return left[0] < right[0] ? -1 : 1;
	});
	toSort.sortIndices = [];
	for (var j = 0; j < toSort.length; j++) 
	{
		toSort.sortIndices.push(toSort[j][1]);
		toSort[j] = toSort[j][0];
	}
	return toSort;
}

function setActiveExclusively(elem)
{
	$('.active').removeClass('active');
	elem.addClass('active');
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

function removeIndexes(arr, rmIdx)
{
	rmIdx.sort(function(a,b){return a - b});
	rmIdx.reverse();
	var narr = arr.slice();

	for (var i = 0; i < rmIdx.length; i++)
	{
		narr.splice(rmIdx[i], 1);
	}

	return narr;
}

function calculateStar(centerX, centerY, arms, outerRadius, innerRadius)
{
   var results = "";

   var angle = Math.PI / arms;

   for (var i = 0; i < 2 * arms; i++)
   {
      // Use outer or inner radius depending on what iteration we are in.
      var r = (i & 1) == 0 ? outerRadius : innerRadius;

      var currX = centerX + Math.cos(i * angle) * r;
      var currY = centerY + Math.sin(i * angle) * r;

      // Our first time we simply append the coordinates, subsequet times
      // we append a ", " to distinguish each coordinate pair.
      if (i == 0)
      {
         results = currX + "," + currY;
      }
      else
      {
         results += ", " + currX + "," + currY;
      }
   }

   return results;
}

function buildConfidenceTable(results)
{
	if (typeof results == 'undefined' || Object.keys(results) == 0)
	{
		return;
	}

	var krakenEnabled = "krakenLabels" in results[Object.keys(results)[0]];

	$('#confidences').append('<tr><th>ID</th><th>Contamination<br/>status</th><th>Sample</th><th>CC</th><th>Dip</th>' + (krakenEnabled ? '<th>Kraken</th>' : '') + '</tr>');

	for (var i in results)
	{
		var kraken = '';
		if (krakenEnabled)
		{
			var krakenNumUnknown = 0;
			$.each(results[i].krakenLabels, function (idx, val)
			{
				if (val === "unknown")
				{
					krakenNumUnknown++;
				}
			});
			var krakenUnique = results[i].krakenLabels.filter(onlyUnique);
			var numSpecies = krakenUnique.length;
			kraken = krakenNumUnknown > 0 ? '&ge; ' + numSpecies : '' + numSpecies;
		}

		var contProbCC = results[i].oneshot.confidenceCC;
		var contProbDip = results[i].oneshot.confidenceDip;
		var status = results[i].oneshot.contaminationStatus;

		$('#confidences').append('<tr>' +
			'<td>' + results[i].id + '</td>' +
			'<td class="' + status + '">&nbsp;</td>' +
			'<td class="selectable dataConf">' + i + '</td>' +
			'<td class="selectable ccConf">conf = ' + contProbCC.toFixed(2) + '</td>' +
			'<td class="selectable dipConf">conf = ' + contProbDip.toFixed(2) + '</td>' +
			(krakenEnabled ? '<td class="selectable kraken"><span class="number">' + kraken + '</span> species</td>' : '') +
			'</tr>');

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

function showData(id, dataMat, labels, tooltips, greyedOut, sixteenS, width, height, padding)
{
	var svg = d3.select("#scatter").html("").append("svg").attr("width", width).attr("height", height);

	var doGrey = typeof greyedOut != 'undefined' && greyedOut.length == labels.length;
	var do16S = typeof sixteenS != 'undefined' && sixteenS.length == labels.length;

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

	var shapes = svg.selectAll("circle").data(dataMat).enter();

	shapes.append("circle")
		  .attr("cx", function(d) { return xScale(d[0]); })
		  .attr("cy", function(d) { return yScale(d[1]); })
		  .attr("r", function(d) {return 3; })
		  .style("fill", function(d, i) {return (doGrey && greyedOut[i]) ? greyColor : colors[labels[i] % colors.length]; })
		  .on('mouseover', tip.show)
		  .on('mouseout', tip.hide);

	if (do16S)
	{
		shapes.append("svg:polygon")
			  .select(function(d, i){ return sixteenS[i] ? this : null; })
			  .attr("id", "star_1")
			  .attr("visibility", "visible")
			  .attr("points", function(d) { return calculateStar(xScale(d[0]), yScale(d[1]), 8, 50, 6); })
			  .style("stroke", highlightColor)
			  .style("stroke-width", 1)
			  .style("fill", highlightColor)
			  .style("fill-opacity", .5)
			  .on('mouseover', tip.show)
			  .on('mouseout', tip.hide)
			  .on('click', function(d, i) { window.location.href = 'export/' + id + '-' + i + '.16s'; });
	}


}

function showVisualization()
{
	var x = results[selectedFasta];

	var clustAnaResult = x.oneshot;

	var dataMat = clustAnaResult[selectedReduction];

	var labels;
	var greyedOut = new Array();
	var sixteenS = new Array();
	var outliers = new Array();
	var additionalInfo = '';

	if (highlight16S && "contains16S" in x && x.contains16S.length > 0)
	{
		sixteenS = x.contains16S;
	}

	if (selectedLabels === 'fasta')
	{
		// labels = numericLabels(x.fastaLabels);
		labels = Array.apply(null, Array(x.fastaLabels.length)).map(Number.prototype.valueOf, -1); // no labels / black color
		additionalInfo = 'Size: ' + (x.stats.numBasepairs/1000000).toFixed(2) + ' Mbp, GC-content: ' + (x.stats.gcContent*100).toFixed(2) + ' %';
	} else if(selectedLabels === 'cc')
	{
		labels = clustAnaResult.clustCC.labels;
		outliers = clustAnaResult.clustCC.outlierClusters;
	} else if(selectedLabels === 'dip')
	{
		if (selectedReduction === 'dataPca')
		{
			if (!selectedNumClusters)
			{
				selectedNumClusters = clustAnaResult.contaminationStatus === 'contaminated' ? clustAnaResult.numClustPca : 1;
				$('#numClusters' + selectedNumClusters).prop('checked', true);
			}
			labels = clustAnaResult.clustsPca[selectedNumClusters-1].labels;
			outliers = clustAnaResult.clustsPca[selectedNumClusters-1].outlierClusters;
		} else if (selectedReduction === 'dataSne')
		{
			if (!selectedNumClusters)
			{
				selectedNumClusters = clustAnaResult.contaminationStatus === 'contaminated' ? clustAnaResult.numClustSne : 1;
				$('#numClusters' + selectedNumClusters).prop('checked', true);
			}
			labels = clustAnaResult.clustsSne[selectedNumClusters-1].labels;
			outliers = clustAnaResult.clustsSne[selectedNumClusters-1].outlierClusters;
		}
	} else if(selectedLabels === 'kraken')
	{
		labels = numericLabels(x.krakenLabels);
		for (var i = 0; i < x.krakenLabels.length; i++)
		{
			if (x.krakenLabels[i] === 'unknown')
			{
				greyedOut.push(true);
			} else
			{
				greyedOut.push(false);
			}
		}
		additionalInfo = 'Bacterial background: ' + (x.krakenBacterialBackground*100).toFixed(2) + '%';
	}

    var tooltips = results[selectedFasta].fastaLabels.slice();
    for (var i = 0; i < tooltips.length; i++)
    {
        var len = x.stats.contigLength[tooltips[i]];
        var gc = x.stats.contigGcContent[tooltips[i]];
        tooltips[i] = tooltips[i] + '<br>Size: ' + (len/1000).toFixed(2) + ' kbp<br>GC-content: ' + (gc*100).toFixed(2) + ' %';
    }
	if(selectedLabels === 'kraken')
	{
		tooltips = x.krakenLabels;
	}


	if (!showOutliers && typeof outliers != 'undefined' && outliers.length > 0)
	{
		var n = dataMat.length;
		var rmIdx = new Array();

		for (var i = 0; i < n; i++)
		{
			if ($.inArray(labels[i], outliers) >= 0)
			{
				rmIdx.push(i);
			}
		}
		rmIdx = arrayUnique(rmIdx);

		dataMat = removeIndexes(dataMat, rmIdx);
		labels = removeIndexes(labels, rmIdx);
		tooltips = removeIndexes(tooltips, rmIdx);
		sixteenS = removeIndexes(sixteenS, rmIdx);
	}

	showData(x.id, dataMat, labels, tooltips, greyedOut, sixteenS, width, height, padding);
	updateExport(labels, greyedOut);

	$('#additionalInfo').html(additionalInfo);
}

function exportClusterFasta(clusterLabel, fasta, selectedLabels, selectedReduction, selectedNumClusters)
{
	var contigNames = {};
	if(selectedLabels !== 'dip')
	{
		contigNames = clusterinfo[fasta][selectedLabels][clusterLabel];
	} else
	{
		contigNames = clusterinfo[fasta][selectedLabels][selectedReduction][selectedNumClusters][clusterLabel];
	}

	// sort by length
	lengths = new Array();
	for (var c in contigNames)
	{
		len = results[fasta].stats.contigLength[contigNames[c]];
		lengths.push(len);	
	}
	sortedLengths = arraySortGetIndices(lengths);
	sortIdx = sortedLengths.sortIndices;
	sortIdx.reverse();
	var sortedContigNames = new Array();
	for (var i in sortIdx)
	{
		sortedContigNames.push(contigNames[sortIdx[i]]);	
	}

	var str = "";
	for (var c in sortedContigNames)
	{
		str += ">" + sortedContigNames[c] + "\n";
		str += inputcontigs[sortedContigNames[c]] + "\n";
	}

	var blob = new Blob([str], {type: "text/plain;charset=utf-8"});
	saveAs(blob, "export.fasta");	
}

function updateExport(labels, greyedOut)
{
	if (selectedLabels === 'fasta')
	{
		$('#export').hide();
		return;
	}

	// find unknown label
	var unknownLbl = -1;
	for (var i in labels)
	{
		if (greyedOut[i])
		{
			unknownLbl = labels[i];
			break;
		}
	}

	var uniqueLabels = arrayUnique(labels).sort();

	$('#exportColors').html('');
	for (var i in uniqueLabels)
	{
		var lbl = uniqueLabels[i];
		var color = (lbl == unknownLbl) ? greyColor : colors[lbl % colors.length];
		var dataParams = "{ 'clusterLabel': '" + lbl + "', 'fasta': '" + results[selectedFasta].fasta + "', 'selectedLabels': '" + selectedLabels + "', 'selectedReduction': '" + selectedReduction + "', 'selectedNumClusters': '" + selectedNumClusters + "' }";
		var exportLink = $('<a data-params="' + dataParams + '" style="color: ' + color + ';" href="#">&#x25CF;</a>').click(function() 
		{
			var me = $(this), data = me.data('params');
			data = data.replace(/\'/g, '"');
			var obj = jQuery.parseJSON(data);
    		exportClusterFasta(obj.clusterLabel, obj.fasta, obj.selectedLabels, obj.selectedReduction, obj.selectedNumClusters);
		});
		$('#exportColors').append('<span></span>').append(exportLink);		
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
