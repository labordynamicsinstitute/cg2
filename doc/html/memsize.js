function memsize() {

var ncells=calc.ncells.value;
var npers=calc.npers.value;
var nfirm=calc.nfirm.value;
var ncov=calc.ncov.value;

if (ncells == '') {alert('Number of cells value must be entered.');
   calc.ncells.focus();return false;}
if (ncells != parseInt(ncells)) {alert('Number of cells value must be numeric');
   calc.ncells.focus();return false;}

if (npers == '') {alert('Number of persons value must be entered.');
   calc.npers.focus();return false;}
if (npers != parseInt(npers)) {alert('Number of persons value must be numeric');
   calc.npers.focus();return false;}

if (nfirm == '') {alert('Number of firms value must be entered.');
   calc.nfirm.focus();return false;}
if (nfirm != parseInt(nfirm)) {alert('Number of firms value must be numeric');
   calc.nfirm.focus();return false;}

if (ncov == '') {alert('Number of covariates value must be entered.');
   calc.ncov.focus();return false;}
if (ncov != parseInt(ncov)) {alert('Number of covariates value must be numeric');
   calc.ncov.focus();return false;}

var ncells=Math.floor(ncells);
var npers=Math.floor(npers);
var nfirm=Math.floor(nfirm);
var ncov=Math.floor(ncov);

var ncoef=ncov+npers+nfirm-1;
var ncoef2=ncov+nfirm-1;
var imem=(4*(2*ncells))/(1024*1024);
var idmem=(8*(ncov+2*npers+nfirm+ncells+2*ncov*ncov+npers*ncov+nfirm*ncov+5*ncoef2 +2*ncoef))/(1024*1024);
var total=imem+idmem;

var total=Math.ceil(total);

calc.total.value=total;

return false;

}
