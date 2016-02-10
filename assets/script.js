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
	});

	$('.ccConf').click(function() {
		selectedFasta = $(this).parent().find('td.dataConf').text();
		selectedLabels = 'cc';
		updateBootStrapsSelect();
		showVisualization();
	});
	
	$('.dipConf').click(function() {
		selectedFasta = $(this).parent().find('td.dataConf').text();
		selectedLabels = 'dip';
		updateBootStrapsSelect();
		showVisualization();
	});

	$('.kraken').click(function() {
		selectedFasta = $(this).parent().find('td.dataConf').text();
		selectedLabels = 'kraken';
		updateBootStrapsSelect();
		showVisualization();
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

	if (!selectedFasta)
	{
		$('#confidences').hide();
		$('#scatterContainer').html('<b>No results to show, yet. Reload page to refresh...<b>');
		return;
	}

	updateBootStrapsSelect();
	showVisualization();

});
