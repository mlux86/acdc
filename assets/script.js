var colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"];

var width = 500;
var height = 500;
var padding = 25;

var selectedFasta = '';
var selectedLabels = 'fasta';
var selectedReduction = 'dataSne';
var selectedData = 'oneshot';

$(document).ready(function() 
{

	buildConfidenceTable(results);

	$('.dataConf').click(function() {
		selectedFasta = $(this).text();
		selectedLabels = 'fasta';
		updateBootStrapsSelect();
		showVisualization();
		setBoldExclusively($(this));
	});

	$('.ccConf').click(function() {
		selectedFasta = $(this).parent().find('td.dataConf').text();
		selectedLabels = 'cc';
		updateBootStrapsSelect();
		showVisualization();
		setBoldExclusively($(this));
	});
	
	$('.dipConf').click(function() {
		selectedFasta = $(this).parent().find('td.dataConf').text();
		selectedLabels = 'dip';
		updateBootStrapsSelect();
		showVisualization();
		setBoldExclusively($(this));
	});

	$('.kraken').click(function() {
		selectedFasta = $(this).parent().find('td.dataConf').text();
		selectedLabels = 'kraken';
		updateBootStrapsSelect();
		showVisualization();
		setBoldExclusively($(this));
	});

	$('.radioReduction').click(function() {
		selectedReduction = $(this).val();
		showVisualization();
	});

	$('#bootstraps').change(function() {
		selectedData = $(this).val();
		showVisualization();
	});

	selectedFasta = $('#confidences').find('td.dataConf:first').text();
	setBoldExclusively($('#confidences').find('td.dataConf:first'));

	if (!selectedFasta)
	{
		$('#confidences').hide();
		$('#scatterContainer').html('<b>No results to show, yet. Reload page to refresh...<b>');
		return;
	}

	updateBootStrapsSelect();	
	showVisualization();

    $(window).scroll(stickyScatter);
    stickyScatter();

});
