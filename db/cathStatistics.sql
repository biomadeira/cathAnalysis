/*create database cathStatistics32;*/
create database if not exists cathStatistics32;
drop table if exists cathStatistics32.cathSuperfamilies;
drop table if exists cathStatistics32.pdbChains;
drop table if exists cathStatistics32.pdbNumberOfChains;

create table cathStatistics32.cathSuperfamilies(cathcode text, domainsInSuperfamily int);
create table cathStatistics32.pdbChains(pdbID text , chain char(1));
create table cathStatistics32. pdbNumberOfChains(pdbID text, chains int);

insert into cathStatistics32.cathSuperfamilies(cathcode, domainsInSuperfamily)  
	select cathcode, count(cathcode) 
	from sflexfit32.CDDFDomains 
	group by cathcode;

insert into cathStatistics32.pdbChains(pdbID, chain)  
select distinct (pdbID), chain from sflexfit32.CDDFDomains;

insert into cathStatistics32.pdbNumberOfChains(pdbID, chains) 
select pdbID, count(chain) 
from cathStatistics32.pdbChains group by  pdbID;
