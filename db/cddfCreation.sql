create database if not exists sflexfit32;
SET @@foreign_key_checks = 0;

drop table if exists   sflexfit32.CDDFDomains;
drop table  if exists  sflexfit32.CDDFSegments;
drop table  if exists  sflexfit32.domainSegments;

create table sflexfit32.CDDFDomains (
	FORMAT    text,
	DOMAIN    char(50) collate latin1_general_cs ,
	pdbID     text,
	chain     char(1),
	VERSION   text,
	VERDATE   date, 
	NAME      text, 
	SOURCE    text,
	CATHCODE  text,
	CLASS     text,
	ARCH      text,
	TOPOL     text,
	HOMOL     text,
	DLENGTH   int,
	DSEQH     text,
	DSEQS     text,
	NSEGMENTS int,
	primary key(DOMAIN)
);

create table sflexfit32.CDDFSegments(
	segment   char(50) collate latin1_general_cs ,
	start     int,
	stop      int, 
	SLENGTH   int,
	SSEQH     text,
	SSEQS     text,
	startInsertion char(1),
	stopInsertion char(2),
	primary key(segment)
);

create table sflexfit32.domainSegments(
	domain  char(50) collate latin1_general_cs ,
	segment char(50) collate latin1_general_cs ,
	foreign key(domain)  references sflexfit32.CDDFDomains(DOMAIN) ,
	foreign key(segment) references sflexfit32.CDDFSegments(segment),
	primary key (domain, segment)
);
SET @@foreign_key_checks = 1;
